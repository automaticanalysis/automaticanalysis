% Eddy corrects image distorsions, head movements using affine registration
% to a reference volume (T2 image)

function [aap resp]=aamod_diffusion_eddycorrect(aap,task,subjind,diffsessind)
global aaworker
resp='';

switch task
    case 'report'
    case 'doit'
        
        % Get nii filenames from stream
        diffinput=aas_getfiles_bystream(aap,'diffusion_session',[subjind diffsessind],'diffusion_data');
        
        % Apply eddy using the reference image (T2 image) which is typycally 0
        % (this is the volume with a b-value of 0)
        [pth nme ext]=aas_fileparts(diffinput);
        eddyoutput=fullfile(pth,[nme '_eddy' ext]);
        
        cmd=sprintf('eddy_correct %s %s 0',diffinput,eddyoutput);
        [s w]=aas_runfslcommand(aap,cmd);
        if s 
            aas_log(aap,true,sprintf('Error executing\n  %s\nof\n%s',cmd,w));
        end  
        
        % Describe outputs
        aap=aas_desc_outputs(aap,'diffusion_session',[subjind diffsessind],'diffusion_data',eddyoutput);   
       
end
end

