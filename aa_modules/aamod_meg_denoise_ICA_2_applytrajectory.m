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
        
        thresholding = 1:size(ICA.ica{1}.temprem,1); % default
        if isfield(aap.tasklist.currenttask.settings,'thresholding') && ~isempty(aap.tasklist.currenttask.settings.thresholding)
            thresholding = aap.tasklist.currenttask.settings.thresholding;
        end
        
        %% Define toremove
        PCA_dim = size(ICA.ica{1}.weights,1);
        toremove = cell(1,numel(ICA.ica));
        for m=1:numel(ICA.ica)
            switch aap.tasklist.currenttask.settings.toremove
                case 'temp' % define artefacts based on temporal correlation only
                    toremove{m} = intersect(1:PCA_dim,unique(cat(2,ICA.ica{m}.temprem{thresholding,:})));
                case 'spat' % define artefacts based on spatial correlation only
                    toremove{m} = intersect(1:PCA_dim,unique(cat(2,ICA.ica{m}.spatrem{thresholding,:})));
                case 'both' % define artefacts based on both
                    toremove{m} = intersect(1:PCA_dim,unique(cat(2,ICA.ica{m}.bothrem{thresholding,:})));
                case 'all' % define artefacts across all references (DEFAULT)
                    toremove{m} = ICA.ica{m}.allrem;
            end
                
            weights = ICA.ica{m}.weights;
            iweights  = pinv(ICA.ica{m}.weights); % weights will be in output file from previous step
            finalics  = setdiff(1:PCA_dim,toremove{m}); % to
            TraMat{m} = iweights(:,finalics) * weights(finalics,:);
        end
        
        %% RUN
        S = []; S.D = D;
        achans = cat(2,ICA.chans{:});
        for c = 1:length(achans)
            S.montage.labelorg{c} = D.chanlabels{achans(c)};
        end
        S.montage.labelnew = S.montage.labelorg;
        S.montage.tra      = blkdiag(TraMat{:});
        S.keepothers       = 1;
        S.prefix           = 'M';
        D = spm_eeg_montage(S);
        
        %% Outputs
        outfname = [S.prefix basename(infname_meg)];
        aap=aas_desc_outputs(aap,subj,sess,'meg',char([outfname '.dat'],[outfname '.mat']));
end