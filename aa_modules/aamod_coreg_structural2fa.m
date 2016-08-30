% AA module - coregister structural to diffusion images (dti_FA)
% Coregistration of structural to dti_FA output by realignment
% Rhodri Cusack BMI, Western, Canada, 2013

function [aap,resp] = aamod_coreg_structural2fa(aap,task,subjInd,diffsess)

resp='';

switch task
	case 'report' 
    case 'doit'
        global defaults;
        flags = defaults.coreg;
        % check local structural directory exists
        structfn=aas_getfiles_bystream(aap,'subject',subjInd,'structural');
        
        dsesspath=aas_getpath_bydomain(aap,'diffusion_session',[subjInd diffsess]);
        
        [pth nme ext]=fileparts(structfn);
        diffsess_structfn=fullfile(dsesspath,[nme ext]);
        copyfile(structfn,diffsess_structfn);
        
        fafn=aas_getfiles_bystream(aap,'diffusion_session',[subjInd diffsess],'dti_FA');

        % Register...
        VG = aas_spm_vol(fafn,true); % ...this
        VF = aas_spm_vol(diffsess_structfn,true); % ...to this
        
        % Replace with unpacked name
        structfn=VF.fname;
        
        % do coregistration
        x  = spm_coreg(VG, VF,flags.estimate);
        
        M  = inv(spm_matrix(x));
          
        spm_get_space(diffsess_structfn, M*spm_get_space(diffsess_structfn));
       
        aap = aas_desc_outputs(aap,'diffusion_session',[subjInd diffsess],'structural', diffsess_structfn);

        % Save graphical output - this will now be done by report task
        try
            figure(spm_figure('FindWin', 'Graphics'));
        catch
            figure(1);
        end
        print('-djpeg','-r75',fullfile(aas_getsubjpath(aap, subjInd),'diagnostic_aamod_coreg_structural2fa'));

	case 'checkrequirements'
        
end
end


