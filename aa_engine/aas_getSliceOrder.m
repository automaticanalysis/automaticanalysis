function aap = aas_getSliceOrder(aap,V, dicomHdr)

if nargin<1
    error('We need an aap structure')
elseif nargin<2
    error('We need a sample volume, as we don''t know the slice number')
elseif nargin<3
    error('We need a dicom header, as we don''t know the slice ordering')
end


if ~isfield(dicomHdr, 'sliceorder')
    str =  dicomHdr.Private_0029_1020;
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
    dicomHdr.sliceorder = sliceorder;

end

switch(dicomHdr.sliceorder)
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
        
        if isfield(dicomHdr, 'Private_0019_1029')
            [~, aap.tasklist.currenttask.settings.sliceorder] = sort(dicomHdr.Private_0019_1029);
            dicomHdr.sliceorder = 'custom (determined from field Private_0019_1029)';
        else
            error('BAD ORDER! Check your slice order, and/or set it manually!');
        end
end

fprintf('\n Your sequence has %d slices in %s order\n', V.dim(3), dicomHdr.sliceorder);

if isfield(aap.tasklist.currenttask.settings, 'slicetime') && isempty(aap.tasklist.currenttask.settings.slicetime)
    if ~isfield(aap.tasklist.currenttask.settings, 'TRs') || ...
            isempty(aap.tasklist.currenttask.settings.TRs)
        % NOTE, this will not work for a 3D sequence
        % Reason 1) A 3D sequence does not have slice order
        % Reason 2) The RepetitionTime is not actually the Volume TR
        aap.tasklist.currenttask.settings.TRs = dicomHdr.RepetitionTime;
    end
    aap.tasklist.currenttask.settings.slicetime=aap.tasklist.currenttask.settings.TRs/V.dim(3);
end
