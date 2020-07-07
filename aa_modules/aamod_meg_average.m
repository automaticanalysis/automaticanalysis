function [aap, resp] = aamod_meg_average(aap,task,varargin)

resp='';

switch task
    case 'report'
        
    case 'doit'
        %% Initialise
        [junk, SPMtool] = aas_cache_get(aap,'spm');
        SPMtool.addCollection('meeg');
        
        switch aap.tasklist.currenttask.domain
            case 'subject'
                localroot = aas_getsubjpath(aap,varargin{1});
            case 'meeg_session'
                localroot = aas_getsesspath(aap,varargin{:});
        end
        infname = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,cell2mat(varargin),'meg'); infname = infname(1,:);        
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
        S.robust = aas_getsetting(aap,'robust');
        S.prefix = 'm';
        D = spm_eeg_average(S);

        %% Outputs
        outfname = fullfile(localroot,[S.prefix basename(infname)]); % specifying output filestem
        fname(D,[outfname '.mat']);
        D.save;

        SPMtool.rmCollection('meeg');

        aap=aas_desc_outputs(aap,aap.tasklist.currenttask.domain,cell2mat(varargin),'meg',char([outfname '.dat'],[outfname '.mat']));        
    case 'checkrequirements'
        if ~aas_cache_get(aap,'spm'), aas_log(aap,true,'SPM is not found'); end
end
end