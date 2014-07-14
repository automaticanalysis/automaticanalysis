% AA module - extended coregistration 
% Coregistration of structural to T1 template
% 1) Reorient base to structural(output 'structural' from aamod_coreg_extended_1)
% 2) Coregister base to structural (using xfm mat called
% 't1totemplate_xfm from aamod_coreg_extended_1)
% 3) Apply all (1 and 2) to inputs.

function [aap,resp]=aamod_coreg_extended_2(aap,task,subj,sess)

resp='';

switch task
%     case 'report' % [TA]
%         if ~exist(fullfile(aas_getsubjpath(aap,subj),['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_structural2meanepi.jpg']),'file')
%             aas_fsl_coreg_diag(aap,subj);
%         end
%         fdiag = dir(fullfile(aas_getsubjpath(aap,subj),'diagnostic_*.jpg'));
%         for d = 1:numel(fdiag)
%             aap = aas_report_add(aap,subj,'<table><tr><td>');
%             aap=aas_report_addimage(aap,subj,fullfile(aas_getsubjpath(aap,subj),fdiag(d).name));
%             aap = aas_report_add(aap,subj,'</td></tr></table>');
%         end
    case 'doit'

        global defaults
        flags = defaults.coreg;
        
        %% Parse inputstreams
        inputs = aap.tasklist.currenttask.inputstreams.stream;
        iEmpty = [];
        for i = 1:numel(inputs)
            if isstruct(inputs{i}), inputs{i} = inputs{i}.CONTENT; end
            if ~aas_stream_has_contents(aap,inputs{i}), iEmpty(end+1) = i; end
        end
        inputs(iEmpty) = [];
        switch aap.tasklist.currenttask.domain
            case 'diffusion_session'
                ibase = cell_index(inputs,'nodif');
                basearg = {aap.tasklist.currenttask.domain,[subj,sess],inputs{ibase}};
            case 'session'
                ibase = cell_index(inputs,'meanepi');
                basearg = {subj,inputs{ibase}};
        end
        inputs(ibase) = [];
        inputs(cell_index(inputs,'structural')) = [];
        inputs(cell_index(inputs,'t1totemplate_xfm')) = [];
        
        %% 0) Check that the templates we need exist!
        % Get the T1 template
        sTimg = fullfile(spm('dir'), aap.directory_conventions.T1template);
        if ~exist(sTimg, 'file')
            aas_log(aap, true, sprintf('Couldn''t find template T1 image %s.', sTimg));
        end
        
        % Get the EPI template
        eTimg = fullfile(spm('dir'), fullfile(fileparts(aap.directory_conventions.T1template),'EPI.nii'));
        if ~exist(eTimg, 'file')
            aas_log(aap, true, sprintf('Couldn''t find template EPI image %s.', eTimg));
        end
        
        % Check local structural directory exists
        
        Simg = aas_getfiles_bystream(aap,subj,'structural');
        if size(Simg,1) > 1
            aas_log(aap, false, sprintf('Found more than 1 structural images, using structural %d', ...
                aap.tasklist.currenttask.settings.structural));
        end

        %% 1) Base to T1 template (reorient)
         
          % Look for base
        
          mEPIimg = aas_getfiles_bystream(aap,basearg{:});
                
            if size(mEPIimg,1) > 1
                aas_log(aap, false, 'Found more than 1 base images, using first.');
                mEPIimg = deblank(mEPIimg(1,:));
            end
         
          % Look for xfm t1totemplate
        
          load(aas_getfiles_bystream(aap,subj,'t1totemplate_xfm'));
  
        % Set the new space for the base
        
           Me = inv(spm_matrix(xfm));
        
           spm_get_space(mEPIimg, Me*spm_get_space(mEPIimg));
          
            fprintf(['\tBase to template realignment parameters:\n' ...
            '\tx: %0.4f   y: %0.4f   z: %0.4f   p: %0.4f   r: %0.4f   j: %0.4f\n'], ...
            xfm(1), xfm(2), xfm(3), xfm(4), xfm(5), xfm(6))
        
        %% 2) Base to Structural (coregister)
            
        % Coregister base to structural
        x = spm_coreg(spm_vol(deblank(Simg(aap.tasklist.currenttask.settings.structural,:))), ...
            spm_vol(mEPIimg(1,:)), ...
            flags.estimate);
        Mf = inv(spm_matrix(x));
        
        % Set the new space for the base
        MM = spm_get_space(mEPIimg);
        spm_get_space(mEPIimg, Mf*MM);
        
        fprintf(['\tmean EPI to structural realignment parameters:\n' ...
            '\tx: %0.4f   y: %0.4f   z: %0.4f   p: %0.4f   r: %0.4f   j: %0.4f\n'], ...
            x(1), x(2), x(3), x(4), x(5), x(6))
            
        %% 3) Now apply this transformation to all the input images
        % The base will already be in the space required for the
        % individual inputs. Hence, we can...
        
        % Again, get space of mean functional
        MM = spm_get_space(mEPIimg(1,:));
        
        EPIimg = cell(size(inputs));
        % Locate all the inputs we want to coregister
        for i = 1:numel(inputs)
            EPIimg{i} = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,[subj,sess],inputs{i});
            fprintf('\nCoregistering %s image(s) for session: %s\n', inputs{i}, aas_getdirectory_bydomain(aap,aap.tasklist.currenttask.domain,sess))
            % For each image, apply the space of the base image
            for e = 1:size(EPIimg{i},1)
                % Apply the space of the coregistered base to the
                % remaining iunputs (safest solution!)
                spm_get_space(deblank(EPIimg{i}(e,:)), MM);
            end
        end
%%  Some Diagnostic Images
%              mriname = aas_prepare_diagnostic(aap,subj);
%          
%              spm_check_registration(strvcat( ...
%              sTimg, ... % Get template T1
%              deblank(Simg(aap.tasklist.currenttask.settings.structural,:)),... % Get structural
%              mEPIimg, ... % Get mean EPI across sessions
%              EPIimg{sess}(1,:))) % Get first image of last session EPI
%          
%              Outline of structural!
%              spm_ov_reorient('context_init', 2)
%          
%              print('-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
%              [mfilename '__' mriname '.jpeg']));
         
%% Diagnostic VIDEO
%          if aap.tasklist.currenttask.settings.diagnostic
%              % Realignment params
%              defs = aap.spm.defaults.realign;
%              
%              % ...flags to pass to routine to create resliced images
%              % (spm_reslice)
%              resFlags = struct(...
%                  'interp', defs.write.interp,...       % interpolation type
%                  'wrap', defs.write.wrap,...           % wrapping info (ignore...)
%                  'mask', defs.write.mask,...           % masking (see spm_reslice)
%                  'which', 1,...     % what images to reslice
%                  'mean', 0);           % write mean image
%              
%              % Get resliced mean EPI
%              [mEPIpth, mEPIfn, mEPIext] = fileparts(deblank(mEPIimg(aap.tasklist.currenttask.settings.structural,:)));
%              spm_reslice(strvcat(Simg, mEPIimg), resFlags);
%              
%              Ydims = {'X', 'Y', 'Z'};
%              for d = 1:length(Ydims)
%                  aas_image_avi( fullfile(mEPIpth, ['r' mEPIfn mEPIext]), ...
%                  Simg, ...
%                  fullfile(aap.acq_details.root, 'diagnostics', [mfilename '__' mriname '_' Ydims{d} '.avi']), ...
%                  d, ... % Axis
%                  [800 600], ...
%                  2); % Rotations
%              end
%              try close(2); catch; end
%              delete(fullfile(mEPIpth, ['r' mEPIfn mEPIext]))
%          end
         
        %% Describe the outputs
        
        aap = aas_desc_outputs(aap,basearg{:},mEPIimg);
        
        for i = 1:numel(inputs)
             aap = aas_desc_outputs(aap,aap.tasklist.currenttask.domain,[subj,sess],inputs{i},EPIimg{i});
        end
    case 'checkrequirements'
        aas_log(aap,0,'No need to trim or skull strip structural\n' );
end