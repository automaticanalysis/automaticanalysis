% AA module
% Runs a premasking of the structural (good to do before using aamod_bet)
% 1) Coregisters the T1 template to the structural
% 2) Reslices this T1 template to the structural
% 3) Masks the structural by the T1 template matrix, cutting out neck

function [aap,resp]=aamod_bet_premask(aap,task,subj)

resp='';

switch task
    case 'domain'
        resp='subject';  % this module needs to be run once per subject
        
    case 'description'
        resp='SPM5 align';
        
    case 'summary'
        subjpath=aas_getsubjpath(subj);
        resp=sprintf('Align %s\n',subjpath);
        
    case 'report'
        
    case 'doit'
        
        Simg = aas_getfiles_bystream(aap,subj,'structural');
        
        % Which file is considered, as determined by the structural parameter!
        if size(Simg,1) > 1
            Simg = deblank(Simg(aap.tasklist.currenttask.settings.structural, :));
            aas_log(aap,false,sprintf('\tWARNING: Several structurals found, considering: %s', Simg))
        end
        
        Sdir = fileparts(Simg);
        
        %% 0) Check that the templates we need exist!
        % Get the template
        % perhaps a sub-dir of SPM?
        sTimg = fullfile(aap.directory_conventions.spmdir,...
            aap.directory_conventions.T1template);
        if ~exist(sTimg,'file')
            % perhaps not  must check this way around because exist produces
            % false positives with relative paths...
            sTimg = aap.directory_conventions.T1template;
        end
        if ~exist(sTimg, 'file')
            aas_log(aap, true, sprintf('Couldn''t find template T1 image %s.', sTimg));
        end
        
        %% 1) Structural to T1 template
        aas_log(aap,false,'Coregistering the T1 template to structural')
        
        global defaults %#ok<TLEV>
        flags = defaults.coreg;
        
        % Copy template to structural location
        copyfile(sTimg, fullfile(Sdir, 'T1.nii'));
        sTimg = fullfile(Sdir, 'T1.nii');
        
        % Coregister template to Structural
        x = spm_coreg(spm_vol(Simg), spm_vol(sTimg), flags.estimate);
        
        % Set entire template matrix to 1
        V = spm_vol(sTimg);
        Y = spm_read_vols(V);
        Y(:) = 1;
        spm_write_vol(V,Y);

        % Set the new space for the template
        MM = spm_get_space(sTimg);
        spm_get_space(sTimg, spm_matrix(x)\MM);
        
        %% 2) Then reslice the Template
        
        aas_log(aap,false,'Reslicing T1 template to structural')
        % Get realignment defaults
        defs = aap.spm.defaults.realign;
        
        % Flags to pass to routine to create resliced images
        % (spm_reslice)
        resFlags = struct(...
            'interp', defs.write.interp,...       % interpolation type
            'wrap', defs.write.wrap,...           % wrapping info (ignore...)
            'mask', defs.write.mask,...           % masking (see spm_reslice)
            'which', 1,...     % what images to reslice
            'mean', 0);           % write mean image
        
        % Reslice
        spm_reslice(strvcat(Simg, sTimg), resFlags);
        
        %% 3) Mask the Structural image using the T1 template
        
        aas_log(aap,false,'Mask structural with resliced T1 template')
        
        M = spm_read_vols(spm_vol(fullfile(Sdir, 'rT1.nii')));
        M = M>0;
        
        % Mask structural
        V = spm_vol(Simg);
        Y = spm_read_vols(V);
        Y = Y.*M;
        spm_write_vol(V,Y);
        
        delete(fullfile(Sdir, 'T1.nii'));
        delete(fullfile(Sdir, 'rT1.nii'));
        
        %% DIAGNOSTIC IMAGE
        % Save graphical output to common diagnostics directory
        if ~exist(fullfile(aap.acq_details.root, 'diagnostics'), 'dir')
            mkdir(fullfile(aap.acq_details.root, 'diagnostics'))
        end

        try
            %% Draw structural image...
            spm_check_registration(Simg)
            
            %% Diagnostic VIDEO of masks
            spm_orthviews('reposition', [0 0 0])
            
            try figure(spm_figure('FindWin', 'Graphics')); catch; figure(1); end;
            set(gcf,'PaperPositionMode','auto')
            print('-djpeg','-r75',fullfile(aap.acq_details.root, 'diagnostics', ...
                [mfilename '__' aap.acq_details.subjects(subj).subjname '.jpeg']));
        catch
        end
        
        %% DESCRIBE OUTPUTS!
        
        % Structural image after BETting
        aap=aas_desc_outputs(aap,subj,'structural', Simg);
        
end
