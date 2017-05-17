% AA module
% Runs MELODIC on all sessions of each single subject
% This automatically transforms the 3D data into 4D data as well
% [NOTE: This function may become obsolete later on]

function [aap,resp]=aamod_melodic(aap,task,subj)

resp='';

switch task
        
    case 'report'
        
    case 'doit'
        
        spaced_EPIimg = [];
        for sess = aap.acq_details.selected_sessions
            % Let us use the native space...
            EPIimg = aas_getfiles_bystream(aap,subj,sess,'epi');
            
            for e = 1:size(EPIimg,1)
                spaced_EPIimg = [spaced_EPIimg EPIimg(e,:) ' '];
            end
        end
        
        %% CONCATENATE THE DATA...
        aas_log(aap,false,'Concatenating the data')
        
        data4D = fullfile(aas_getsubjpath(aap,subj), sprintf('4Ddata_%s.nii', aap.acq_details.subjects(subj).subjname));
        
        [junk, w]=aas_runfslcommand(aap, ...
            sprintf('fslmerge -t %s %s', ...
            data4D, ...
            spaced_EPIimg));
        
        %% RUN MELODIC
        aas_log(aap,false,'Running MELODIC')
        
        outDir = fullfile(aas_getsubjpath(aap,subj), 'MELODIC');
        if ~exist(outDir, 'dir')
            mkdir(outDir)
        end
        
        [junk, w]=aas_runfslcommand(aap, ...
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
                if ~D(d).isdir && isempty(strfind(D(d).name(1), '.'))
                    melodicFiles = strvcat(melodicFiles, fullfile(fldrCurr, D(d).name));
                else
                    % It is one of the . or .. folders
                end
            end
        end
        
        aap=aas_desc_outputs(aap,subj,'melodic', melodicFiles);
        
end