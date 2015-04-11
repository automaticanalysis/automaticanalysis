% AA module - extended coregistration 
% Coregistration of structural to T1 template


function [aap,resp]=aamod_coreg_extended_1(aap,task,subj)

resp='';

switch task
    case 'doit'

        global defaults
        flags = defaults.coreg;
        
        %% 0) Check that the templates we need exist!
        % Get the T1 template
        sTimg = aap.directory_conventions.T1template;
        if (~exist(sTimg,'file')) % try in SPM
            if sTimg(1) ~= '/', sTimg = fullfile(fileparts(which('spm')),sTimg); end
        end
        if ~exist(sTimg,'file'),
            aas_log(aap, true, sprintf('Couldn''t find template T1 image %s.', sTimg));
        end  
               
        %% 1) Structural to T1 template
        % Check local structural directory exists
        
        Simg = aas_getfiles_bystream(aap,subj,'structural');
        if size(Simg,1) > 1
            aas_log(aap, false, sprintf('Found more than 1 structural images, using structural %d', ...
                aap.tasklist.currenttask.settings.structural));
        end
        
        % Coregister T1 to T1 template
        xfm = spm_coreg(spm_vol(sTimg), ...
            spm_vol(deblank(Simg(aap.tasklist.currenttask.settings.structural,:))), ...
            flags.estimate);
        fn = fullfile(aas_getsubjpath(aap,subj),'t1totemplate_xfm.mat');
        save (fn,'xfm')
        Ms = inv(spm_matrix(xfm));
        
        % Set the new space for the structural
        for d = 1:size(Simg,1)
            MM = spm_get_space(deblank(Simg(d,:)));
            spm_get_space(deblank(Simg(d,:)), Ms*MM);
        end
        
        fprintf(['\tstructural to template realignment parameters:\n' ...
            '\tx: %0.4f   y: %0.4f   z: %0.4f   p: %0.4f   r: %0.4f   j: %0.4f\n'], ...
            xfm(1), xfm(2), xfm(3), xfm(4), xfm(5), xfm(6))
        
%         
        %% Describe the outputs
        
        aap = aas_desc_outputs(aap,subj,'structural',Simg);
        aap = aas_desc_outputs(aap, subj,'t1totemplate_xfm',fn);
        
      
end