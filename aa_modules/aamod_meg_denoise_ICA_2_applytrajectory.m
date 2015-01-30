function [aap, resp] = aamod_meg_denoise_ICA_2_applytrajectory(aap,task,subj,sess)

resp='';

switch task
    case 'report'
        
    case 'doit'
        %% Initialise
        infname_meg = aas_getfiles_bystream(aap,'meg_session',[subj sess],'meg'); infname_meg = infname_meg(1,:);
        infname_ica = aas_getfiles_bystream(aap,'meg_session',[subj sess],'meg_ica');
        
        
        D=spm_eeg_load(infname_meg);  % INPUTSTREAM from Convert
        ICA = load(infname_ica);
        
        %% RUN
        S = []; S.D = D;
        achans = cat(2,ICA.chans{:});
        for c = 1:length(achans)
            S.montage.labelorg{c} = D.chanlabels{achans(c)};
        end
        S.montage.labelnew = S.montage.labelorg;
        S.montage.tra      = blkdiag(ICA.TraMat{:});
        S.keepothers       = 1;
        S.prefix           = 'M';
        D = spm_eeg_montage(S);
        
        %% Outputs
        outfname = [S.prefix basename(infname_meg)];
        aap=aas_desc_outputs(aap,subj,sess,'meg',char([outfname '.dat'],[outfname '.mat']));
end