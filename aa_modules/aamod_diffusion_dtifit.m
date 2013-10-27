% DTIFIT (FSL software) fits a diffusion tensor model at each voxel.
%Note that dtifit is not neccessary in order to run  probabilistic
%tractrography (which depends on the output of BEDPOSTX)

function [aap resp]=aamod_diffusion_dtifit(aap,task,subjind,diffsessind)
resp='';

switch task
    case 'report'
    case 'doit'
        
        % Get nii filenames from stream
        diffinput=aas_getfiles_bystream(aap,'diffusion_session',[subjind diffsessind],'diffusion_data');
        bvals=aas_getfiles_bystream(aap,'diffusion_session',[subjind diffsessind],'bvals');
        bvecs=aas_getfiles_bystream(aap,'diffusion_session',[subjind diffsessind],'bvecs');
        betmask=aas_getfiles_bystream(aap,'diffusion_session',[subjind diffsessind],'BETmask');
        
        % Find which line of betmask contains the brain mask
        for betind=1:size(betmask,1)
            if strfind(betmask(betind,:),'bet_nodif_brain_mask')
                break
            end;
        end;        
        betmask=betmask(betind,:);
        
        % Get diffusion data streams and apply dtifit

        fslext=aas_getfslext(aap);
        dsesspth= aas_getpath_bydomain(aap,'diffusion_session',[subjind,diffsessind]);
        out=fullfile(dsesspth,'dti');
        cmd=sprintf('dtifit -b %s -r %s  -k %s  -m %s -o %s',bvals,bvecs,diffinput,betmask,out)
        [s w]=aas_runfslcommand(aap,cmd);
        if (s)
            aas_log(aap,true,sprintf('Error executing\n  %s\nof\n%s',cmd,w));
        end;
        
        % Now describe outputs
        outstreams=aap.tasklist.currenttask.outputstreams;        
        for outind=1:length(outstreams.stream)
            aap=aas_desc_outputs(aap,'diffusion_session',[subjind,diffsessind],outstreams.stream{outind},[outstreams.stream{outind} fslext]);
        end;
        
end
end



