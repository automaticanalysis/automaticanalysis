% AA module - slice timing
% [aap,resp]=aamod_slicetiming(aap,task,subj,sess)
% Corrects the slice time difference in 2D EPI sequences
% Do not use on 3D EPI (it will not work!)
% Rhodri Cusack MRC CBU 2004 based on original by Matthew Brett
% Batching slice timing in SPM5
% Tibor Auer MRC CBU Cambridge 2012-2013

function [aap,resp]=aamod_slicetiming(aap,task,subj,sess)
resp='';

switch task
    case 'domain'
        resp='session';   % this module needs to be run once per session
    case 'description'
        resp='SPM8 slice timing';
    case 'summary'
        resp=sprintf('Perform slice timing, with TR %f and time to acquire single slice %f\n',aap.tasklist.currenttask.settings.TRs, aap.tasklist.currenttask.settings.slicetime);
    case 'report'
    case 'doit'
        
        % get the subdirectories in the main directory
        dirn = aas_getsesspath(aap,subj,sess);
        
        % get files in this directory
        imgs=aas_getimages_bystream(aap,subj,sess,'epi');
        
        % get information from first file
        first_img = deblank(imgs(1,:));
        DICOMHEADERS = load(aas_getimages_bystream(aap,subj,sess,'epi_dicom_header'));
        hdr = DICOMHEADERS.DICOMHEADERS{1};
        V = spm_vol(first_img);
        if isfield(aap.options, 'NIFTI4D') && aap.options.NIFTI4D % 4D support [TA]
            V = V(1);
        end
        
        % retrieve stuff from DICOM header
        if aap.tasklist.currenttask.settings.autodetectSO == 1
            % Get the headers from the file, so that we don't have to guess...            
            if isnumeric(hdr.sliceorder) && ~isempty(hdr.slicetimes) % exact info
                % save for outputs
                sliceorder=hdr.sliceorder;
                refslice=aap.tasklist.currenttask.settings.refslice;
                
                aap.tasklist.currenttask.settings.sliceorder = hdr.slicetimes*1000;
                aap.tasklist.currenttask.settings.refslice = aap.tasklist.currenttask.settings.sliceorder(refslice);
                sl_times = 0;
            else
                if ~isempty(hdr.sliceorder) && isfield(hdr,'Private_0029_1020')
                    aap = aas_getSliceOrder(aap, V, hdr);            
                end
            end
        end
        if isempty(aap.tasklist.currenttask.settings.TRs)
            aap.tasklist.currenttask.settings.TRs=hdr.RepetitionTime/1000;
        end
        
        % Sets slice time information
        if exist('sl_times','var') % exact slicetiming info
            % value 1 is 0
            % value 2 is TR
            sl_times = [0 aap.tasklist.currenttask.settings.TRs];
        else
            % value 1 is time to acquire one slice
            % value 2 is time between beginning of last slice and beginning of first slice of next volume
            if isempty(aap.tasklist.currenttask.settings.slicetime) % rough estimate only
                aap.tasklist.currenttask.settings.slicetime=aap.tasklist.currenttask.settings.TRs/V.dim(3);
            end
            if (max(aap.tasklist.currenttask.settings.sliceorder)>V.dim(3))
                aas_log(aap,1,'aap.tasklist.currenttask.settings.sliceorder seems to contain values higher than the number of slices!\n');
            end
            sl_times = [aap.tasklist.currenttask.settings.slicetime aap.tasklist.currenttask.settings.slicetime+(aap.tasklist.currenttask.settings.TRs-aap.tasklist.currenttask.settings.slicetime*V.dim(3))];
            % outputs
            sliceorder=aap.tasklist.currenttask.settings.sliceorder;
            refslice=aap.tasklist.currenttask.settings.refslice;
        end
        
        % do slice timing correction
        spm_slice_timing(imgs,aap.tasklist.currenttask.settings.sliceorder,aap.tasklist.currenttask.settings.refslice ,sl_times);
        
        % Describe outputs
        rimgs=[];
        for k=1:size(imgs,1);
            [pth nme ext]=fileparts(imgs(k,:));
            rimgs=strvcat(rimgs,['a' nme ext]);
        end
        aap = aas_desc_outputs(aap,subj,sess,'epi',rimgs);
        
        sliceorderfn=fullfile(dirn,'sliceorder.mat');
        save(sliceorderfn,'sliceorder','refslice');
        aap = aas_desc_outputs(aap,subj,sess,'sliceorder',sliceorderfn);
        
    case 'checkrequirements'
        if (length(aap.tasklist.currenttask.settings.sliceorder)==0) && aap.tasklist.currenttask.settings.autodetectSO == 0
            aas_log(aap,1,'To avoid catastrophe, slice order no longer takes on a default value, and must be specified manually in user script.\nFor descending sequential 32 slices add aap.tasksettings.aamod_slicetiming.sliceorder=[32:-1:1];\n');
        end
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
return;
end

function aap = aas_getSliceOrder(aap,V, dicomHdr)

if nargin<1
    aas_log(aap,true,'We need an aap structure')
elseif nargin<2
    aas_log(aap,true,'We need a sample volume, as we don''t know the slice number')
elseif nargin<3
    aas_log(aap,true,'We need a dicom header, as we don''t know the slice ordering')
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
        aas_log(aap,false,'WARNING: Ensure your interleaved order is correct!');
    otherwise
        
        if isfield(dicomHdr, 'Private_0019_1029')
            [junk, aap.tasklist.currenttask.settings.sliceorder] = sort(dicomHdr.Private_0019_1029);
            dicomHdr.sliceorder = 'custom (determined from field Private_0019_1029)';
        else
            aas_log(aap,true,'BAD ORDER! Check your slice order, and/or set it manually!');
        end
end

aas_log(aap,false,sprintf('INFO: Your sequence has %d slices in %s order', V.dim(3), dicomHdr.sliceorder));

if isfield(aap.tasklist.currenttask.settings, 'slicetime') && isempty(aap.tasklist.currenttask.settings.slicetime)
    if ~isfield(aap.tasklist.currenttask.settings, 'TRs') || ...
            isempty(aap.tasklist.currenttask.settings.TRs)
        % NOTE, this will not work for a 3D sequence
        % Reason 1) A 3D sequence does not have slice order
        % Reason 2) The RepetitionTime is not actually the Volume TR
        % NOTE, this will not work for sparse imaging
        aap.tasklist.currenttask.settings.TRs = dicomHdr.RepetitionTime;
    end
    aap.tasklist.currenttask.settings.slicetime=aap.tasklist.currenttask.settings.TRs/V.dim(3);
end
end