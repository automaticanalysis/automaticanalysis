function [aap, resp] = aamod_meg_epochs(aap,task,subj,sess)

resp='';

switch task
    case 'report'
        
    case 'doit'
        %% Initialise
        [junk, SPMtool] = aas_cache_get(aap,'spm');
        SPMtool.addCollection('meeg');
        
        sessdir = aas_getsesspath(aap,subj,sess);
        infname = aas_getfiles_bystream(aap,'meeg_session',[subj sess],'meeg'); infname = infname(1,:);
        S.D = spm_eeg_load(infname);
        
        S.bc = aap.tasklist.currenttask.settings.baselinecorrection;
        S.timewin = aap.tasklist.currenttask.settings.timewindow;
        
        %% Conditions
        subjmatches=strcmp(aap.acq_details.subjects(subj).subjname,{aap.tasklist.currenttask.settings.condition.subject});
        sessmatches=strcmp(aap.acq_details.meeg_sessions(sess).name,{aap.tasklist.currenttask.settings.condition.session});
        % If no exact spec found, try session wildcard, then subject
        % wildcard, then wildcard for both
        if (~any(sessmatches & subjmatches))
            sesswild=strcmp('*',{aap.tasklist.currenttask.settings.condition.session});
            if (any(sesswild & subjmatches))
                sessmatches=sesswild;
            else
                subjwild=strcmp('*',{aap.tasklist.currenttask.settings.condition.subject});
                if (any(sessmatches & subjwild))
                    subjmatches=subjwild;
                else
                    subjmatches=subjwild;
                    sessmatches=sesswild;
                end
            end
        end
        
        % Should now have just one model spec
        conditionnum = sessmatches & subjmatches;        
        
        if ~isempty([aap.tasklist.currenttask.settings.condition(conditionnum).event.trlshift]) % data-specified
            S.trialdef = aap.tasklist.currenttask.settings.condition(conditionnum).event;
        else  % user-specified
            S.trl = aap.tasklist.currenttask.settings.condition(conditionnum).event.eventvalue;
            if ~isempty(aap.tasklist.currenttask.settings.condition(conditionnum).event.conditionlabel)
                S.conditionlabels = aap.tasklist.currenttask.settings.condition(conditionnum).event.conditionlabel;
            end
        end
           
        %% Run
        S.prefix = 'e';
        D = spm_eeg_epochs(S);

        %% Outputs
        outfname = fullfile(sessdir,[S.prefix basename(infname)]); % specifying output filestem
        fname(D,[outfname '.mat']);
        D.save;
        
        SPMtool.rmCollection('meeg');
        
        aap=aas_desc_outputs(aap,subj,sess,'meeg',char([outfname '.dat'],[outfname '.mat']));        
    case 'checkrequirements'
        if ~aas_cache_get(aap,'spm'), aas_log(aap,true,'SPM is not found'); end
end
end