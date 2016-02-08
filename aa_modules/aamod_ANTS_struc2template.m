% AA module
% Use the Anatomical Transformation Toolbox to normalise the structural to
% a template image

function [aap,resp]=aamod_ANTS_struc2template(aap,task,subj)

resp='';

switch task
    case 'doit'
        
        warning off
        
        %% Template image
        % [AVG] Changed to allow specification of any T1 template, does not
        % need to be in the SPM folder any more...
        sTimg = aap.directory_conventions.T1template;
        if ~exist(sTimg, 'file')
            aas_log(aap, true, sprintf('Couldn''t find template T1 image %s.', Timg));
        end
        
        %% Get structural
        % [AVG] Modified the way we get the structural, to be more aa4-like
        Simg = aas_getfiles_bystream(aap,subj,'structural');
        % Cheap and cheerful way of ensuring only one file is considered!
        if size(Simg,1) > 1
            Simg = deblank(Simg(aap.tasklist.currenttask.settings.structural,:));
            aas_log(aap,0,sprintf('Found more than one structural so using:\n%s',Simg));
        end
        % Get structural directory for this subject
        [Spth, Sfn, Sext] = fileparts(Simg);
        
        %% Use ANTS to normalise them!
        
        % Set the ANTS path
        setenv('ANTSPATH', aap.directory_conventions.ANTSdir)
        ANTSpath = [fullfile(getenv('ANTSPATH'), 'bin', 'ANTS') ' '];
        warpANTSpath = [fullfile(getenv('ANTSPATH'), 'bin', 'WarpImageMultiTransform') ' '];
        
        % What we get out...
        outfiles = '-o ants.nii ';
        
        % Set up maxiterations
        maxiterations = '';
        for m = 1:length(aap.tasklist.currenttask.settings.maxiterations)
            maxiterations = [maxiterations num2str(aap.tasklist.currenttask.settings.maxiterations(m)) 'x'];
        end
        maxiterations = ['-i ' maxiterations(1:end-1) ' '];
        
        % Dimension number (always 3 for structural)
        Ndim = [num2str(3) ' '];
        
        % Regularisation...
        if ~isempty(aap.tasklist.currenttask.settings.regularisation)
            regularisation = aap.tasklist.currenttask.settings.regularisation;
        else
            regularisation = '';
        end
        
        % SyN transformation...
        if ~isempty(aap.tasklist.currenttask.settings.SyN)
            SyN = ['-t SyN[' num2str(aap.tasklist.currenttask.settings.SyN) '] '];
        else
            SyN = '';
        end
        
        metrics = '';
        for m = 1:9
            if isfield(aap.tasklist.currenttask.settings, ['metric' num2str(m)])
                tmpM = aap.tasklist.currenttask.settings.(['metric' num2str(m)]);
                tmpW = num2str(aap.tasklist.currenttask.settings.(['weight' num2str(m)]));
                tmpP = aap.tasklist.currenttask.settings.(['parameters' num2str(m)]);
                if isnumeric(tmpP)
                    tmpP = num2str(tmpP);
                end
                
                metrics = [ metrics ...
                    '-m ' tmpM '[' sTimg ',' Simg ',' tmpW ',' tmpP '] '];
            else
                break
            end
        end
        
        % Any extra options?...
        if ~isempty(aap.tasklist.currenttask.settings.SyN)
            extraoptions = aap.tasklist.currenttask.settings.extraoptions;
        else
            extraoptions = '';
        end
        
        ANTS_command = [ ANTSpath Ndim outfiles maxiterations SyN metrics extraoptions];
        
        cd(Spth)
        
        % Run ANTS
        fprintf('Running ANTS using command:\n')
        fprintf([ANTS_command '\n'])
        [s w] = aas_shell(ANTS_command);
        
        warpANTS_command = [ warpANTSpath Ndim ... % dimension number
            Simg ' ' fullfile(Spth, ['w' Sfn Sext]) ... % moving image & output
            ' -R ' sTimg ' '... % reference image
            fullfile(Spth, 'antsWarp.nii')]; % transform
        if exist(fullfile(Spth,'antsAffine.txt'), 'file')
            warpANTS_command = [warpANTS_command ' antsAffine.txt']; % and affine, if this exists...
        end    
        
        [s w] = aas_shell(warpANTS_command);
        
        %% Describe outputs
        aap=aas_desc_outputs(aap,subj,'structural', fullfile(Spth,['w' Sfn Sext]));
        outANTS = strvcat( ...
        fullfile(Spth,'antsWarp.nii'), ...
        fullfile(Spth, 'antsInverseWarp.nii'));
        if exist(fullfile(Spth,'antsAffine.txt'), 'file')
            outANTS = strvcat(outANTS, fullfile(Spth,'antsAffine.txt'));
        end
        aap=aas_desc_outputs(aap,subj,'ANTs', outANTS);
        
        % Diagnostic image?
        % Save graphical output to common diagnostics directory
        if ~exist(fullfile(aap.acq_details.root, 'diagnostics'), 'dir')
            mkdir(fullfile(aap.acq_details.root, 'diagnostics'))
        end
        
        %% Draw native template
        spm_check_registration(strvcat( ...
            fullfile(Spth,['w' Sfn Sext]), ...
            sTimg))
        
        %% Diagnostic VIDEO of segmentations
        aas_checkreg_avi(aap, subj, 2)
        
        spm_orthviews('reposition', [0 0 0])
        
        try figure(spm_figure('FindWin', 'Graphics')); catch; figure(1); end;
        print('-djpeg','-r75',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' aap.acq_details.subjects(subj).subjname '.jpeg']));
        
end
end