% AA module - extended coregistration of EPI to structural
% Coregistration of structural to mean EPI output by realignment in 3 steps
% 1) Coregister Structural to T1 template
% 2) Coregister mean EPI to EPI template
% 3) Coregister mean EPI to Structural
% 4) Apply transformation matrix of mean EPI to all EPIs

function [aap,resp]=aamod_coreg_extended(aap,task,subj)

resp='';

switch task
    case 'report' % [TA]
        if ~exist(fullfile(aas_getsubjpath(aap,subj),['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_structural2meanepi.jpg']),'file')
%             aas_fsl_coreg_diag(aap,subj);
        end
        fdiag = dir(fullfile(aas_getsubjpath(aap,subj),'diagnostic_*.jpg'));
        for d = 1:numel(fdiag)
            aap = aas_report_add(aap,subj,'<table><tr><td>');
            aap=aas_report_addimage(aap,subj,fullfile(aas_getsubjpath(aap,subj),fdiag(d).name));
            aap = aas_report_add(aap,subj,'</td></tr></table>');
        end
    case 'doit'

        global defaults
        flags = defaults.coreg;
        
        
        %% 0) Check that the tamplates we need exist!
        % Get the template
        sTimg = fullfile(spm('dir'), aap.directory_conventions.T1template);
        if ~exist(sTimg, 'file')
            aas_log(aap, true, sprintf('Couldn''t find template T1 image %s.', sTimg));
        end
        
        % Get the template
        eTimg = fullfile(fileparts(sTimg), 'EPI.nii');
        if ~exist(eTimg, 'file')
            aas_log(aap, true, sprintf('Couldn''t find template EPI image %s.', eTimg));
        end
        
        %% 1) Structural to T1 template
        % Check local structural directory exists
        
        Simg = aas_getfiles_bystream(aap,subj,'structural');
        if size(Simg,1) > 1
            aas_log(aap, false, sprintf('Found more than 1 structural images, using structural %d', ...
                aap.tasklist.currenttask.settings.structural));
        end
        
        % Coregister T1 to template
        x = spm_coreg(spm_vol(sTimg), ...
            spm_vol(deblank(Simg(aap.tasklist.currenttask.settings.structural,:))), ...
            flags.estimate);
        Ms = inv(spm_matrix(x));
        
        % Set the new space for the structural
        for d = 1:size(Simg,1)
            MM = spm_get_space(deblank(Simg(d,:)));
            spm_get_space(deblank(Simg(d,:)), Ms*MM);
        end
        
        fprintf(['\tstructural to template realignment parameters:\n' ...
            '\tx: %0.4f   y: %0.4f   z: %0.4f   p: %0.4f   r: %0.4f   j: %0.4f\n'], ...
            x(1), x(2), x(3), x(4), x(5), x(6))
        
        %% 2) Mean Functional to EPI template
        
        % Look for mean functional
        mEPIimg = aas_getfiles_bystream(aap,subj,'meanepi');
                
        if size(mEPIimg,1) > 1
            aas_log(aap, false, 'Found more than 1 mean functional images, using first.');
            mEPIimg = deblank(mEPIimg(1,:));
        end
        
        % Coregister mean functional to template
        x = spm_coreg(spm_vol(eTimg), spm_vol(mEPIimg), flags.estimate);
        Me = inv(spm_matrix(x));
        
        % Set the new space for the mean functional
        MM = spm_get_space(mEPIimg(1,:));
        spm_get_space(mEPIimg, Me*MM);
        
        fprintf(['\tmean EPI to template realignment parameters:\n' ...
            '\tx: %0.4f   y: %0.4f   z: %0.4f   p: %0.4f   r: %0.4f   j: %0.4f\n'], ...
            x(1), x(2), x(3), x(4), x(5), x(6))
        
        %% 3) Mean Functional to Structural
            
        % Coregister mean EPI to structural
        x = spm_coreg(spm_vol(deblank(Simg(aap.tasklist.currenttask.settings.structural,:))), ...
            spm_vol(mEPIimg(1,:)), ...
            flags.estimate);
        Mf = inv(spm_matrix(x));
        
        % Set the new space for the mean EPI
        MM = spm_get_space(mEPIimg);
        spm_get_space(mEPIimg, Mf*MM);
        
        fprintf(['\tmean EPI to structural realignment parameters:\n' ...
            '\tx: %0.4f   y: %0.4f   z: %0.4f   p: %0.4f   r: %0.4f   j: %0.4f\n'], ...
            x(1), x(2), x(3), x(4), x(5), x(6))
        
        %% 4) Now apply this transformation to all the EPI images
        % The mean EPI will already be in the space required for the
        % individual EPIs. Hence, we can...
        
        % Again, get space of mean functional
        MM = spm_get_space(mEPIimg(1,:));
        
        EPIimg = cell(size(aap.acq_details.sessions));
        % Locate all the EPIs we want to coregister
        for sess = aap.acq_details.selected_sessions
            EPIimg{sess} = aas_getfiles_bystream(aap,subj,sess,'epi');
            
            % For each image, apply the space of the mean EPI image
            fprintf('\nCoregistering images for session: %s\n', aas_getsessdesc(aap,subj,sess))
            for e = 1:size(EPIimg{sess},1)
                % Apply the space of the coregistered mean EPI to the
                % remaining EPIs (safest solution!)
                spm_get_space(deblank(EPIimg{sess}(e,:)), MM);
            end
        end
        
        %% Some diagnostic images
        subjname = aas_prepare_diagnostic(aap,subj);
        
        spm_check_registration(strvcat( ...
            sTimg, ... % Get template T1
            deblank(Simg(aap.tasklist.currenttask.settings.structural,:)),... % Get structural
            mEPIimg, ... % Get mean EPI across sessions
            [EPIimg{sess}(1,:) ',1'])) % Get first image of last session EPI
        
        % Outline of structural!
        spm_ov_reorient('context_init', 2)
        
        f = spm_figure('GetWin','Graphics');
        set(f,'Renderer','zbuffer');
        print(f,'-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' subjname '.jpeg']));
        
        %% Diagnostic VIDEO
        if aap.tasklist.currenttask.settings.diagnostic
            % Realignment params
            defs = aap.spm.defaults.realign;
            
            % ...flags to pass to routine to create resliced images
            % (spm_reslice)
            resFlags = struct(...
                'interp', defs.write.interp,...       % interpolation type
                'wrap', defs.write.wrap,...           % wrapping info (ignore...)
                'mask', defs.write.mask,...           % masking (see spm_reslice)
                'which', 1,...     % what images to reslice
                'mean', 0);           % write mean image
            
            % Get resliced mean EPI
            [mEPIpth, mEPIfn, mEPIext] = fileparts(deblank(mEPIimg(aap.tasklist.currenttask.settings.structural,:)));
            spm_reslice(strvcat(Simg, mEPIimg), resFlags);
            
            Ydims = {'X', 'Y', 'Z'};
            for d = 1:length(Ydims)
                aas_image_avi( fullfile(mEPIpth, ['r' mEPIfn mEPIext]), ...
                Simg, ...
                fullfile(aap.acq_details.root, 'diagnostics', [mfilename '__' subjname '_' Ydims{d} '.avi']), ...
                d, ... % Axis
                [800 600], ...
                2); % Rotations
            end
            try close(2); catch; end
            delete(fullfile(mEPIpth, ['r' mEPIfn mEPIext]))
        end
        
        %% Describe the outputs
        
        aap = aas_desc_outputs(aap,subj,'structural',Simg);
        aap = aas_desc_outputs(aap,subj,'meanepi',mEPIimg);
        
        for sess = aap.acq_details.selected_sessions
            aap = aas_desc_outputs(aap,subj,sess,'epi',EPIimg{sess});
        end
        
    case 'checkrequirements'
        aas_log(aap,0,'No need to trim or skull strip structural\n' );
end