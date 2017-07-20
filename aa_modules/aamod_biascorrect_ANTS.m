% AA module
% Use the Anatomical Transformation Toolbox to normalise the structural to
% a template image

function [aap,resp]=aamod_biascorrect_ANTS(aap,task,subj)

resp='';

switch task
    case 'doit'
        
        %% Get image to bias correct
        % Find out what stream we should use
        inputstream = aap.tasklist.currenttask.inputstreams.stream;
        % And the names of the output streams
        outputstream = aap.tasklist.currenttask.outputstreams.stream;
        
        % Let us get the image we want to bias correct...
        Simg = aas_getfiles_bystream(aap,subj,inputstream{:});
        
        %% Use ANTS to bias correct them
        
        % Set the ANTS path
        setenv('ANTSPATH', aap.directory_conventions.ANTSdir)
        ANTSpath = [fullfile(getenv('ANTSPATH'), 'bin', 'N4BiasFieldCorrection') ' '];
        
        % Dimension number (always 3 for structural)
        Ndim = [num2str(3) ' '];
        
        % Any extra options?...
        if ~isempty(aap.tasklist.currenttask.settings.extraoptions)
            extraoptions = aap.tasklist.currenttask.settings.extraoptions;
        else
            extraoptions = '';
        end
        
        outstream = ''; biasstream = '';
        for d = aap.tasklist.currenttask.settings.structural
            % Get structural directory for this subject
            [Spth, Sfn, Sext] = fileparts(deblank(Simg(d,:)));
            
            outstream = strvcat(outstream, fullfile(Spth, ['m' Sfn Sext]));
            biasstream = strvcat(biasstream, fullfile(Spth, ['bias' Sfn Sext]));
            
            infiles = sprintf('-i %s ', ...
                fullfile(Spth, [Sfn Sext]));
            
            outfiles = sprintf('-o [%s,%s] ', ...
                outstream(d==aap.tasklist.currenttask.settings.structural, :), ...
                biasstream(d==aap.tasklist.currenttask.settings.structural, :));
            
            ANTS_command = [ ANTSpath Ndim infiles outfiles extraoptions];
            
            cd(Spth)
            
            % Run ANTS
            [s w] = aas_shell(ANTS_command);
        end
        
        %% Describe outputs
        aap=aas_desc_outputs(aap,subj,outputstream{:}, outstream);
        
        %% DIAGNOSTIC
        subjname = aas_prepare_diagnostic(aap,subj);
        
        for d = aap.tasklist.currenttask.settings.structural
            %% Draw native template
            spm_check_registration(strvcat( ...
                deblank(Simg(d,:)), ...
                outstream(d==aap.tasklist.currenttask.settings.structural, :), ...
                biasstream(d==aap.tasklist.currenttask.settings.structural, :)))
            
            %% Diagnostic VIDEO of segmentations
            aas_checkreg_avi(aap, subj, 2)
            
            spm_orthviews('reposition', [0 0 0])
            
            print('-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
                [mfilename '__' subjname '_' num2str(d) '.jpeg']));
        end
end
end