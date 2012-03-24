function aap = aas_getSliceOrder(aap, p, s, V)

if nargin<1
    error('We need an aap structure')
elseif nargin<2
    error('We need a participant number')
elseif nargin<3
    error('We need a session number')
elseif nargin<4
    error('We need a sample volume, as we don''t know the slice number')
end


dicomPath = (fullfile(aap.directory_conventions.rawdatadir, ...
    aap.acq_details.subjects(p).mriname));

dicomFold = dir(dicomPath);

dicomFold = dicomFold(aap.acq_details.subjects(p).seriesnumbers(s) + 2).name;

dicomDir = dir(fullfile(dicomPath, dicomFold));

dicomName = fullfile(dicomPath, dicomFold, dicomDir(3).name);

infoD = dicominfo(dicomName);
str = infoD.Private_0029_1020;
xstr = char(str');
n = findstr(xstr, 'sSliceArray.ucMode');
[t, r] = strtok(xstr(n:n+100), '=');
ucmode = strtok(strtok(r, '='));
switch(ucmode)
    case '0x1'
        sliceorder = 'Ascending';
    case '0x2'
        sliceorder = 'Descending';
    case '0x4'
        sliceorder = 'Interleaved';
    otherwise
        sliceorder = 'Order undetermined';
end

switch(sliceorder)
    case 'Ascending'
        [aap.tasklist.currenttask.settings.sliceorder] = 1:1:V.dim(3);
    case 'Descending'
        [aap.tasklist.currenttask.settings.sliceorder] = V.dim(3):-1:1;
    case 'Interleaved'
        % Interleaved order depends on whether slice number is odd or even!
        if mod(V.dim(3),2)
            [aap.tasklist.currenttask.settings.sliceorder] = [1:2:V.dim(3) 2:2:V.dim(3)];
        else
            [aap.tasklist.currenttask.settings.sliceorder] = [2:2:V.dim(3) 1:2:V.dim(3)];
        end
        warning('CAREFUL! Ensure your interleaved order is correct!')
    otherwise
        error('BAD ORDER! Check your slice order, and/or set it manually!')
end

fprintf('\n Your sequence has %d slices in %s order', V.dim(3), sliceorder);

if isempty(aap.tasklist.currenttask.settings.slicetime)
    if ~isfield(aap.tasklist.currenttask.settings, 'TRs') || ...
            isempty(aap.tasklist.currenttask.settings.TRs)
        % NOTE, this will not work for a 3D sequence
        % Reason 1) A 3D sequence does not have slice order
        % Reason 2) The RepetitionTime is not actually the Volume TR
        aap.tasklist.currenttask.settings.TRs = infoD.RepetitionTime;
    end
    aap.tasklist.currenttask.settings.slicetime=aap.tasklist.currenttask.settings.TRs/V.dim(3);
end

