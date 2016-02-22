% AA module - normalisation using normalise or two pass procedure with segment
% [aap,resp]=aamod_norm_noss(aap,task,subj)
% Depending on aap.tasksettings.aamod_norm_noss.usesegmentnotnormalise
% If 0 - classic Normalisation
% If 1 - use Segment with two pass procedure: a first pass to correct
%  for intensity variation across our structural images (probably due to
%  inhomogeneous SNR across space of 12 channel coil); and a second pass
%  to then do the segmentation
% _noss version does not use skull stripping
% subj=subject num
% Rhodri Cusack & Daniel Mitchell MRC CBU 2006
% based on originals by Rik Henson, Matthew Brett

function [aap,resp]=aamod_biascorrect(aap,task,subj)
resp='';

switch task
    case 'report'
        
    case 'doit'
        
        defs =aap.spm.defaults.normalise;
        defs.estimate.weight = '';        
        
        %% Get image to bias correct
        % Find out what stream we should use
        inputstream = aap.tasklist.currenttask.inputstreams.stream;
        % And the names of the output streams
        outputstream = aap.tasklist.currenttask.outputstreams.stream;
        
        % Let us get the image we want to bias correct...
        Simg = aas_getfiles_bystream(aap,subj,inputstream{:});
        
        % Get structural directory for this subject
        [Spth, Sfn, Sext] = fileparts(Simg(aap.tasklist.currenttask.settings.structural(1),:));
        
        mimgfn = '';
        for d = aap.tasklist.currenttask.settings.structural
            %% Set up normalisation, etc.
            % 2 stage process, as proposed by RH, to increased robustness [djm 13/03/06]
            
            %%%%%%%% 1st pass:
            fprintf('Running first pass of norm_noss (get bias corrected structural)\n')
            estopts.regtype='';    % turn off affine:
            out = spm_preproc(Simg(d,:), estopts);
            [sn,isn]   = spm_prep2sn(out);
            
            % only write out attenuation corrected image
            writeopts.biascor = 1;
            writeopts.GM  = [0 0 0];
            writeopts.WM  = [0 0 0];
            writeopts.CSF = [0 0 0];
            writeopts.cleanup = [0];
            spm_preproc_write(sn, writeopts);
            
            % get the filename for the bias-corrected image (which has 'm' prepended)
            [pth, nm, ext] = fileparts(Simg(d,:));
            mimgfn = strvcat(mimgfn, fullfile(pth, sprintf('m%s%s', nm, ext)));
            
            %% Diagnostic image?
            subjname = aas_prepare_diagnostic(aap,subj);
            
            % Draw image pre/post bias correction
            spm_check_registration(strvcat(Simg(d,:), ...
                mimgfn(aap.tasklist.currenttask.settings.structural==d,:)))
            
            spm_orthviews('reposition', [0 0 0])
            
            f = spm_figure('GetWin','Graphics');
            set(f,'Renderer','zbuffer');
            print(f, '-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
                [mfilename '__' subjname '_' num2str(d) '.jpeg']));
        end
        
        %% DESCRIBE OUTPUTS
        aap=aas_desc_outputs(aap,subj, outputstream{:}, mimgfn);

        %% Diagnostic image?
        subjname = aas_prepare_diagnostic(aap,subj);
        
        % Draw image pre/post bias correction
        spm_check_registration(strvcat(Simg, ...
            fullfile(Spth,['m' Sfn Sext])))
        
        spm_orthviews('reposition', [0 0 0])
        
        print('-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' subjname '.jpeg']));
        
end
