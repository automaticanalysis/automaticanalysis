function [aap, resp] = aamod_meg_average(aap,task,subj,sess)

resp='';

switch task
    case 'report'
        
    case 'doit'
        %% Initialise
        sessdir = aas_getsesspath(aap,subj,sess);
        infname = aas_getfiles_bystream(aap,'meg_session',[subj sess],'meg'); infname = infname(1,:);
        D = spm_eeg_load(infname);
        
        % Reorder conditions
        if ~isempty(aap.tasklist.currenttask.settings.conditionorder)
            conditionorder = textscan(aap.tasklist.currenttask.settings.conditionorder,'%s','delimiter',':'); conditionorder = conditionorder{1}';
            if any(~strcmp(sort(D.condlist), sort(conditionorder)))
                c = D.condlist;
                aas_log(aap,true,sprintf('ERROR: Conditions specified here are not in concordance with conditions specified in aamod_meg_epochs: %s',sprintf('\n\t%s',c{:})));
            end
            D = condlist(D, conditionorder);
            D.save
        end
       
        %% Run
        S.D = D;
        S.robust = aap.tasklist.currenttask.settings.robust;
        S.prefix = 'm';
        D = spm_eeg_average(S);

        %% Outputs
        outfname = fullfile(sessdir,[S.prefix basename(infname)]); % specifying output filestem
        fname(D,[outfname '.mat']);
        D.save;
        aap=aas_desc_outputs(aap,subj,sess,'meg',char([outfname '.dat'],[outfname '.mat']));        
end