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
            % value 2 is timeRepetitionTime between beginning of last slice and beginning of first slice of next volume
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
