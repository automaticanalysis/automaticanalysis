function [aap, resp] = aamod_meg_merge(aap,task,subj)

resp='';

switch task
    case 'doit'
        fname = {};
        for sess = aap.acq_details.selected_sessions
            fname{end+1} = aas_getfiles_bystream(aap,'meg_session',[subj sess],'meg');
        end
        
        cd(aas_getsubjpath(aap,subj));
        
        S.D = char(fname);
        S.recode = aas_getsetting(aap,'recode');
        S.prefix = aas_getsetting(aap,'prefix');
        Dout = spm_eeg_merge(S);
        
        %% Outputs
        aap=aas_desc_outputs(aap,'subject',subj,'meg',char(spm_file(Dout.fnamedat,'ext','mat'),Dout.fnamedat));
end