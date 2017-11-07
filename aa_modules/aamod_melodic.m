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
        
        %% Brain mask (optional)
        if aas_stream_has_contents(aap,subj,'firstlevel_brainmask')
            aap.tasklist.currenttask.settings.MELODICoptions = ['-m ' aas_getfiles_bystream(aap,'subject',subj,'firstlevel_brainmask') ' '...
                aap.tasklist.currenttask.settings.MELODICoptions];
        end
        
        %% RUN MELODIC
        aas_log(aap,false,'Running MELODIC')
        
        outDir = fullfile(aas_getsubjpath(aap,subj), 'MELODIC');
        aas_makedir(aap,outDir);
        
        [junk, w]=aas_runfslcommand(aap, ...
            sprintf('melodic -i %s %s -o %s', ...
            data4D, ...
            aap.tasklist.currenttask.settings.MELODICoptions, ...
            outDir));
        
        % Delete 4D file once we finish!
        delete(data4D);
        
        %% DESCRIBE OUTPUTS!
        melodicFiles = spm_select('FPListRec',outDir);
        aap=aas_desc_outputs(aap,subj,'melodic', melodicFiles);
        
end