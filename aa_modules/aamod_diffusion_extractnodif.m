% Extract the reference(s) image(s) (T2 image with b-value of 0), called
% nodif

function [aap resp]=aamod_diffusion_extractnodif(aap,task,subjind,diffsessind)
resp='';

switch task
    case 'report'
    case 'doit'
        
        % Get nii filenames from stream
        difffn=aas_getfiles_bystream(aap,'diffusion_session',[subjind diffsessind],'diffusion_data');
        if size(difffn,1)>1
            aas_log(aap,true,sprintf('expecting a single 4d file but got %d files',length(difffn)));
        end;
        
        % Load up this file
        [V Y]=aas_spm_vol(difffn);
        
        % Pick out the ones with no diffusion
        bvals=spm_load(aas_getfiles_bystream(aap,'diffusion_session',[subjind diffsessind],'bvals'));
        Y=Y(:,:,:,bvals==0);

        % And average them (if there are more than one)
        Y=mean(Y,4);
        
        % Save the output
        Vout=V(1);
        Vout.fname=fullfile(aas_getpath_bydomain(aap,'diffusion_session',[subjind diffsessind]),'nodif.nii');
        spm_write_vol(Vout,Y);
        
        % Describe outputs
        aap=aas_desc_outputs(aap,'diffusion_session',[subjind diffsessind],'nodif',Vout.fname);   
       
end
end

