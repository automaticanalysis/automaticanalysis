function [aap, resp] = aamod_meg_convert(aap,task,subj,sess)

resp='';

switch task
    case 'report'
        
    case 'doit'
        %% Initialise
        sessdir = aas_getsesspath(aap,subj,sess);
        inputstreams = aas_getstreams(aap,'input');
        infname = aas_getfiles_bystream(aap,'meg_session',[subj sess],inputstreams{2});        
        outfname = fullfile(sessdir,['spm12_' basename(infname)]); % specifying output filestem
        
        chan=load(aas_getfiles_bystream(aap,'channellabels'));
        
        S=[];
        S.dataset       = infname;
        S.outfile       = outfname;
        S.save          = 0;
        S.reviewtrials  = 0;
        S.channels      = chan.label;  % RH Could just pass grads + EOG/ECG here (no mags nor STI101)
        S.continuous    = 1;
        S.checkboundary = 0;
        
        D = spm_eeg_convert(S);
        
        % Relevant to triggers, and no triggers at rest!
        tc = find(strcmp(D.chanlabels,'STI101'));
        if ~isempty(tc)
            D(tc,:) = D(tc,:)/1e6;  % 1e6 new scaling introduced by SPM12!
            D = units(D,tc,'V');
            D.save;
        end
        
        %% Outputs
        aap=aas_desc_outputs(aap,subj,sess,'meg',char([outfname '.dat'],[outfname '.mat']));        
end