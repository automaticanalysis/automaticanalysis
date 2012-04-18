% AA module
% Runs MELODIC on all sessions of each single subject
% This automatically transforms the 3D data into 4D data as well
% [NOTE: This function may become obsolete later on]

function [aap,resp]=aamod_melodic(aap,task,p)

resp='';

switch task
        
    case 'report'
        
    case 'doit'
        
        spaced_EPIimg = [];
        for s = aap.acq_details.selected_sessions
            % Let us use the native space...
            EPIimg = aas_getfiles_bystream(aap,p,s,'epi');
            
            for e = 1:size(EPIimg,1)
                spaced_EPIimg = [spaced_EPIimg EPIimg(e,:) ' '];
            end
        end
        
        mriname = strtok(aap.acq_details.subjects(p).mriname, '/');
        
        %% CONCATENATE THE DATA...
        fprintf('\nConcatenating the data')
        
        data4D = fullfile(aas_getsubjpath(aap,p), sprintf('4Ddata_%s.nii', mriname));
        
        [~, w]=aas_runfslcommand(aap, ...
            sprintf('fslmerge -t %s %s', ...
            data4D, ...
            spaced_EPIimg));
        
        %% RUN MELODIC
        fprintf('\nRunning MELODIC')
        
        outDir = fullfile(aas_getsubjpath(aap,p), 'MELODIC');
        if ~exist(outDir, 'dir')
            mkdir(outDir)
        end
        
        [~, w]=aas_runfslcommand(aap, ...
            sprintf('melodic -i %s %s -o %s', ...
            data4D, ...
            aap.tasklist.currenttask.settings.MELODICoptions, ...
            outDir));
        
        % Delete 4D file once we finish!
        unix(['rm ' data4D])
        
        %% DESCRIBE OUTPUTS!
        
        % MAKE A SEPARATE FUNCTION OF THIS SOMETIME?
        melodicFiles = [];
        fldrDir = genpath(outDir);
        % Then recurse inside each directory until you run out of paths
        while ~isempty(strtok(fldrDir, ':'))
            % Get each of the directories made by gendir
            [fldrCurr fldrDir] = strtok(fldrDir, ':');
            % Check it's not a .svn folder
            D = dir(fldrCurr);
            for d = 1:length(D)
                if ~D(d).isdir || isempty(strfind(D(d).name(1), '.'))
                    melodicFiles = strvcat(melodicFiles, fullfile(outDir, D(d).name));
                else
                    % It is one of the . or .. folders
                end
            end
        end
        
        aap=aas_desc_outputs(aap,p,'melodic', melodicFiles);
        
end