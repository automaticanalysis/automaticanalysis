% Extract the reference(s) image(s) (T2 image with b-value of 0), called
% nodif

function [aap resp]=aamod_diffusion_extractnodif(aap,task,subjind,diffsessind)
resp='';

switch task
    case 'report'
    case 'doit'
        
       
        if aas_stream_has_contents(aap,'diffusion_session',[subjind diffsessind])
            % Just a single phase encode direction
            difffn_list={aas_getfiles_bystream(aap,'diffusion_session',[subjind diffsessind],'diffusion_data')};
            bvalfn_list={aas_getfiles_bystream(aap,'diffusion_session',[subjind diffsessind],'bvals')};
        else
            % Top-up style, multiple phase encode directions
            for peind=1:2
                difffn_list{peind}=aas_getfiles_bystream(aap,'diffusion_session_phaseencode_direction',[subjind diffsessind peind],'diffusion_data');
                bvalfn_list{peind}=aas_getfiles_bystream(aap,'diffusion_session_phaseencode_direction',[subjind diffsessind peind],'bvals');
            end;
        end;
        
        % Load up one or more phase encoding directions
        Ytot=[];
        for peind=1:length(difffn_list)
            if size(difffn_list{peind},1)>1
                aas_log(aap,true,sprintf('expecting a single 4d file but got %d files',length(difffn)));
            end;
            % Load up this file
            [V Y]=aas_spm_vol(difffn_list{peind});
            
            % Pick out the ones with no diffusion
            bvals=spm_load(bvalfn_list{peind});
            
            % Select only images with b=0 and average them
            Y=Y(:,:,:,bvals==0);
            Y=mean(Y,4);
            
            if length(difffn_list)>1
                 % Save the separate phase encodes
                Vout=V(1);
                Vout.fname=fullfile(aas_getpath_bydomain(aap,'diffusion_session_phaseencode_direction',[subjind diffsessind peind]),'nodif_phaseencode_direction.nii');
                spm_write_vol(Vout,Y)

                % Describe outputs
                aap=aas_desc_outputs(aap,'diffusion_session_phaseencode_direction',[subjind diffsessind peind],'nodif',Vout.fname);
                
                Ytot(:,:,:,peind)=Y;
            else
                Ytot=Y;
            end;
            
        end;
        
        % Mean across phase encode directions
        Y=mean(Y,4);
        
        % Save the output
        Vout=V(1);
        Vout.fname=fullfile(aas_getpath_bydomain(aap,'diffusion_session',[subjind diffsessind]),'nodif.nii');
        spm_write_vol(Vout,Y);
        
        % Describe outputs
        aap=aas_desc_outputs(aap,'diffusion_session',[subjind diffsessind],'nodif',Vout.fname);
        
end
end

