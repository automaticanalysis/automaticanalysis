% AA module - slice timing
% [aap,resp]=aamod_slicetiming(aap,task,i,j)
% Corrects the slice time difference in 2D EPI sequences
% Do not use on 3D EPI (it will not work!)
% Rhodri Cusack MRC CBU 2004 based on original by Matthew Brett
% Batching slice timing in SPM5

function [aap,resp]=aamod_slicetiming(aap,task,i,j)
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
        dirn = aas_getsesspath(aap,i,j);
        
        % get files in this directory
        % Old style, by prefix (still supported for now)
        %imgs=aas_getimages(aap,i,j,aap.tasklist.currenttask.epiprefix,0,inf);
        % New style, by stream
        imgs=aas_getimages_bystream(aap,i,j,'epi');
        
        if aap.tasklist.currenttask.settings.autodetectSO == 1
            V = spm_vol(deblank(imgs(1,:)));
            aap = aas_getSliceOrder(aap, i, j, V);
        end
        
        % get information from first file
        first_img = deblank(imgs(1,:));
        V = spm_vol(first_img);
        
        % retrieve stuff from DICOM header
        if (length(aap.tasklist.currenttask.settings.TRs)==0)
            DICOMHEADERS=load(fullfile(dirn,'dicom_headers'));
            aap.tasklist.currenttask.settings.TRs=DICOMHEADERS.DICOMHEADERS{1}.RepetitionTime/1000;
        end
        if (length(aap.tasklist.currenttask.settings.slicetime)==0)
            aap.tasklist.currenttask.settings.slicetime=aap.tasklist.currenttask.settings.TRs/V.dim(3);
        end
        % Sets slice time information
        % value 1 is time to acquire one slice
        % value 2 is time between beginning of last slice
        % and beginning of first slice of next volume
        
        if (max(aap.tasklist.currenttask.settings.sliceorder)>V.dim(3))
            aas_log(aap,1,'aap.tasklist.currenttask.settings.sliceorder seems to contain values higher than the number of slices!\n');
        end
        
        sl_times = [aap.tasklist.currenttask.settings.slicetime aap.tasklist.currenttask.settings.slicetime+(aap.tasklist.currenttask.settings.TRs-aap.tasklist.currenttask.settings.slicetime*V.dim(3))];
        
        % do slice timing correction, added refslice [de 200606]
        spm_slice_timing(imgs,aap.tasklist.currenttask.settings.sliceorder,aap.tasklist.currenttask.settings.refslice ,sl_times);
        
        % Describe outputs
        rimgs=[];
        for k=1:size(imgs,1);
            [pth nme ext]=fileparts(imgs(k,:));
            rimgs=strvcat(rimgs,['a' nme ext]);
        end
        sessdir=aas_getsesspath(aap,i,j);
        aap = aas_desc_outputs(aap,i,j,'epi',rimgs);
        
        sliceorder=aap.tasklist.currenttask.settings.sliceorder;
        sliceorderfn=fullfile(dirn,'sliceorder.mat');
        refslice=aap.tasklist.currenttask.settings.refslice;
        save(sliceorderfn,'sliceorder','refslice');
        aap = aas_desc_outputs(aap,i,j,'sliceorder',sliceorderfn);
        
    case 'checkrequirements'
        if (length(aap.tasklist.currenttask.settings.sliceorder)==0) && aap.tasklist.currenttask.settings.autodetectSO == 0
            aas_log(aap,1,'To avoid catastrophe, slice order no longer takes on a default value, and must be specified manually in user script.\nFor descending sequential 32 slices add aap.tasksettings.aamod_slicetiming.sliceorder=[32:-1:1];\n');
        end
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
return;









