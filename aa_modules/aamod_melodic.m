% AA module
% Runs MELODIC on all sessions of each single subject
% This automatically transforms the 3D data into 4D data as well
% [NOTE: This function may become obsolete later on]

function [aap,resp]=aamod_melodic(aap,task,subj)

resp='';

switch task
        
    case 'report'
        html = fullfile(aas_getsubjpath(aap,subj),'MELODIC','report','00index.html');
        
        aap = aas_report_add(aap,subj,'<table><tr>');
        aap = aas_report_add(aap,subj,sprintf('<td><a href="%s">FSL MELODIC report: %s</a><br></td>',html,html));
        aap = aas_report_add(aap,subj,'</tr></table>');
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
        aas_runfslcommand(aap, ...
            sprintf('fslmerge -t %s %s', data4D, spaced_EPIimg));
        
        %% MELODIC options
        opts = aap.tasklist.currenttask.settings.MELODICoptions;
        
        % Brain mask (optional)
        if aas_stream_has_contents(aap,subj,'firstlevel_brainmask')
            opts = [opts ' -m ' aas_getfiles_bystream(aap,'subject',subj,'firstlevel_brainmask')];
        end
        
        % thresholding (optional)
        mmthresh = aas_getsetting(aap,'mmthresh');
        if isempty(mmthresh), mmthresh = 0.5; end % backward compatibility: 0.5 is default
        if mmthresh > 0
            opts = sprintf('%s --mmthresh %1.3f',opts, mmthresh);
        else
            opts = [opts ' --no_mm'];
        end
        
        %% RUN MELODIC
        aas_log(aap,false,'Running MELODIC')
        
        outDir = fullfile(aas_getsubjpath(aap,subj), 'MELODIC');
        aas_makedir(aap,outDir);
        
        aas_runfslcommand(aap, ...
            sprintf('melodic -i %s %s -o %s --report', data4D, opts, outDir));
        
        % Delete 4D file once we finish!
        delete(data4D);
        
        %% DESCRIBE OUTPUTS!
        melodicFiles = spm_select('FPList',outDir);
        aap=aas_desc_outputs(aap,subj,'melodic', melodicFiles);        
end