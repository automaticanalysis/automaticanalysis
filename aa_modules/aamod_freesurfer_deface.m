% Automated Defacing Tools (Freesurfer software)
% It defaces structural (T1) and produces a mask which can be applied to
% any coregistered image

function [aap resp]=aamod_freesurfer_deface(aap,task,subjind)
resp='';

switch task
    case 'report'
    case 'doit'
        % Check templates
        global aa
        if ~isa(aa,'aaClass'), aa = aaClass; end
        
        tmpdir = fullfile(aa.Path,'external','Templates');
        tmp_face = fullfile(tmpdir,'freesurfer_daface_face.gca');
        tmp_skull = fullfile(tmpdir,'freesurfer_daface_talairach_mixed_with_skull.gca');
        if ~exist(tmp_face,'file') || ~exist(tmp_skull,'file')
            aas_log(aap,true,sprintf('Templates required: %s, %s\n',tmp_face, tmp_skull));
        end
        
        % Get input
        sdir = pwd;
        Simg = aas_getfiles_bystream(aap,subjind,'structural'); 
        [p, fn, ext] = fileparts(Simg);
        out = fullfile(p,['defaced_' fn ext]);
        
        % Apply mri_deface
        cd(p);
        cmd = sprintf('mri_deface %s %s %s %s',Simg,tmp_skull,tmp_face,out);
        [s, w]=aas_runFScommand(aap,cmd);
        if (s)
            aas_log(aap,true,sprintf('Error executing\n  %s\nof\n%s',cmd,w));
        end;
        cd(sdir);
        
        % Create mask (for other images)
        inf = spm_vol(out);
        Y = spm_read_vols(inf);
        Y = Y>0;
        out_mask = fullfile(p,['defaced_mask_' fn ext]);
        nifti_write(out_mask,Y,'Defaced mask',inf)
        
        % Now describe outputs
        aap=aas_desc_outputs(aap,subjind,'defaced_structural',out);
        aap=aas_desc_outputs(aap,subjind,'defaced_mask',out_mask);
end
end



