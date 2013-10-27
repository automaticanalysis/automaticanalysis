% AA module - coregister structural to diffusion images (dti_FA)
% Coregistration of structural to dti_FA output by realignment
% Does not require skull stripping any more
% Modified for sparse imaging since prefix for mean is different
% i=subject num
% Rhodri Cusack MRC CBU 2004-6 based on original by Matthew Brett
% 
% Major changes Aug 2010: removed support for central store of structrual
% images. This code was very long in tooth, and unloved.
%
% Tibor Auer MRC CBU Cambridge 2012-2013

function [aap,resp] = aamod_coreg_structural2fa(aap,task,subjInd,diffsess)

resp='';

switch task
	case 'report' 
    case 'doit'
        global defaults;
        flags = defaults.coreg;
        % check local structural directory exists
        structfn=aas_getfiles_bystream(aap,'subject',subjInd,'structural');
        fafn=aas_getfiles_bystream(aap,'diffusion_session',[subjInd diffsess],'dti_FA');

        % Register...
        VG = aas_spm_vol(fafn,true); % ...this
        VF = aas_spm_vol(structfn,true); % ...to this
        
        % Replace with unpacked name
        structfn=VF.fname;
        
        % do coregistration
        x  = spm_coreg(VG, VF,flags.estimate);
        
        M  = inv(spm_matrix(x));
          
        spm_get_space(structfn, M*spm_get_space(structfn));
       
        aap = aas_desc_outputs(aap,'subject',subjInd,'structural', structfn);

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


