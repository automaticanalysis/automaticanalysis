% AA module
% Runs BET (FSL Brain Sextration Toolbox) on structural
% [For best functionality, it is recommended you run this after
% realigSfnnt and before writing the normalised EPI image
% If you do it before estimating the normalisation, make sure you normalise
% to a scull-stripped template, if at all possible!]

function [aap,resp]=aamod_bet_freesurfer(aap,task,subj)

resp='';

switch task
    case 'report'
        
    case 'doit'
        
        % Let us use the native space...
        Simg = aas_getfiles_bystream(aap,subj,'structural');
        
        % Which file is considered, as determined by the structural parameter!
        if size(Simg,1) > 1
            Simg = deblank(Simg(aap.tasklist.currenttask.settings.structural, :));
            for t = 1:length(aap.tasklist.currenttask.settings.structural)
                fprintf('\t%s\n', Simg(t,:))
            end
        end
        
        % Image that we will be using for BET...
        cSimg = deblank(Simg(1,:));
        [Spth Sfn Sext]=fileparts(cSimg);
        
        % Structural image (after BETting, if we mask...)
        bSimg=fullfile(Spth,['bet_' Sfn Sext]);
        
        % Set the Freesurfer path (and to command)
        setenv('FREESURFER_HOME', aap.directory_conventions.freesurferdir)
        
        FSpath = [fullfile(getenv('FREESURFER_HOME'), 'bin', 'mri_watershed') ' '];
        
        FSoptions = aap.tasklist.currenttask.settings.extraoptions;
        
        FScommand = [FSpath FSoptions ' ' cSimg ' ' bSimg];
        
        % Run MRI_Watershed
        fprintf('MRI Watershed\n')
        cd(fileparts(bSimg));
        [s, w] = aas_shell(FScommand,~aap.tasklist.currenttask.settings.verbose);
        
        if aap.tasklist.currenttask.settings.gcut
            % Gcut
            FSpath = [fullfile(getenv('FREESURFER_HOME'), 'bin', 'mri_gcut') ' '];
            FScommand = [FSpath ' ' cSimg ' ' fullfile(Spth, 'gcut.nii')];
            
            % Run MRI_Gcut
            fprintf('MRI Gcut\n')
            cd(fileparts(bSimg));
            [s, w] = aas_shell(FScommand,~aap.tasklist.currenttask.settings.verbose);
        end
        
        %% Make the BET BRAIN MASK
        V = spm_vol(bSimg);
        M = spm_read_vols(V);
        % Mask out non-brain
        M = M > 0;
        fprintf('Watershed mask contains %d voxels\n', sum(M(:)))
        
        if aap.tasklist.currenttask.settings.gcut
            gV = spm_vol(fullfile(Spth, 'gcut.nii'));
            gM = spm_read_vols(gV);
            % Mask out non-brain
            gM = gM > 0;
            cM = and(M, gM);
            delete(fullfile(Spth, 'gcut.nii'));
            fprintf('Graph-cut mask contains %d voxels\n', sum(gM(:)))
            fprintf('Combined mask contains %d voxels\n', sum(cM(:)))
            
            % Test the overlap (Jaccart) between masks...
            if (sum(and(M(:), gM(:))) ./ sum(or(M(:),gM(:)))) > aap.tasklist.currenttask.settings.jaccart
                M = cM;
            elseif aap.tasklist.currenttask.settings.dominant == 0
                % Watershed is the more liberal mask...
                fprintf('\tCombination fails, so we use only Watershed\n')
            elseif aap.tasklist.currenttask.settings.dominant == 1
                % But Graph cut seems to work better with MP2RAGE...
                fprintf('\tCombination fails, so we use only Graph Cut\n')
                M = gM;
            end
        end
        
        % Then write out actual BET mask
        V.fname = fullfile(Spth, ['bet_' Sfn '_brain_mask' Sext]);
        spm_write_vol(V,M);
        % And add the mask to the list...
        outMask = V.fname;
        
        %% MASK the brain(s)
        fprintf('Masking the brain(s) with Brain Mask \n')
        outStruct = '';
        for t = 1:length(aap.tasklist.currenttask.settings.structural)
            % Mask structural
            V = spm_vol(deblank(Simg(t,:)));
            Y = spm_read_vols(V);
            % Mask brain
            Y = Y.*M;
            % Write brain
            [pth nme ext]=fileparts(deblank(Simg(t,:)));
            V.fname = fullfile(pth,['bet_' nme ext]);
            spm_write_vol(V,Y);
            % Add to output...
            outStruct = strvcat(outStruct, V.fname);
        end
        
        %% DESCRIBE OUTPUTS!
        if aap.tasklist.currenttask.settings.maskBrain
            aap=aas_desc_outputs(aap,subj,'structural',outStruct);
        else
           aap=aas_desc_outputs(aap,subj,'structural',Simg); 
        end
        aap=aas_desc_outputs(aap,subj,'BETmask',outMask);
        
        %% DIAGNOSTIC IMAGE
        subjname = aas_prepare_diagnostic(aap,subj);
        
        %% Draw structural image...
        spm_check_registration(Simg)
        
        % Colour the brain Sextracted bit pink
        spm_orthviews('addcolouredimage',1,outStruct(1,:), [0.9 0.4 0.4])
        
        spm_orthviews('reposition', [0 0 0])
        
        print('-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' subjname '.jpeg']));
        
        %% Diagnostic VIDEO of masks
        if aap.tasklist.currenttask.settings.diagnostic
            
            Ydims = {'X', 'Y', 'Z'};
            for d = 1:length(Ydims)
                aas_image_avi(cSimg, ...
                    fullfile(Spth, ['bet_' Sfn '_brain_mask' Sext]), ...
                    fullfile(aap.acq_details.root, 'diagnostics', [mfilename '__' subjname '_' Ydims{d} '.avi']), ...
                    d, ... % Axis
                    [800 600], ...
                    2); % Rotations
            end
            try close(2); catch; end
        end
        
        % Clean up...
        if ~aap.tasklist.currenttask.settings.maskBrain
            delete(outStruct);
        end
end
