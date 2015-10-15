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
        if isempty(intersect(thresholding,[1:4]))
            aas_log(aap,true,'ERROR: thresholding option must be one or more of [1:4]');
        end
        
        %% Define toremove
    
        PCA_dim = size(ICA.ica{1}.weights,1);
        toremove = cell(1,numel(ICA.ica));
        for m=1:numel(ICA.ica)
            toremove{m} = [1:PCA_dim];
            switch aap.tasklist.currenttask.settings.toremove
                case 'temp' % define artefacts based on temporal correlation thresholds only
                    for t = 1:length(thresholding)
                        toremove{m} = intersect(toremove{m},unique(cat(2,ICA.ica{m}.temprem{thresholding(t),:})));
                    end
                case 'spat' % define artefacts based on spatial correlation thresholds only
                    for t = 1:length(thresholding)
                        toremove{m} = intersect(toremove{m},unique(cat(2,ICA.ica{m}.spatrem{thresholding(t),:})));
                    end
                case 'either' % define artefacts based on either spatial OR temporal correlation thresholds
                    tempremove = [1:PCA_dim];
                    for t = 1:length(thresholding)
                        tempremove = intersect(tempremove,unique(cat(2,ICA.ica{m}.temprem{thresholding(t),:})));
                    end
                    spatremove = [1:PCA_dim];
                    for t = 1:length(thresholding)
                        spatremove = intersect(spatremove,unique(cat(2,ICA.ica{m}.spatrem{thresholding(t),:})));
                    end
                    toremove{m} = unique([tempremove spatremove]);
                case 'both' % define artefacts that meet both temporal AND spatial correlation thresholds (DEFAULT)
                % !!!NOTE: This assumes temporal and spatial references were paired!!!
                    if size(ICA.ica{m}.temprem,2) ~= size(ICA.ica{m}.spatrem,2)
                        aas_log(aap,true,'ERROR: Must have same number of temporal and spatial references to use "both" option');
                    else 
                        for r = 1:size(ICA.ica{m}.temprem,2)
                            for t = 1:length(thresholding)                                
                                both = intersect(ICA.ica{m}.temprem{thresholding(t),r},ICA.ica{m}.spatrem{thresholding(t),r});
                                toremove{m} = intersect(toremove{m},both);
                            end
                        end
                    end
            end
            
            toremove{m} = intersect(toremove{m},ICA.ica{m}.varenough); % Additional criterion of sufficient variance (normally ignored, ie varenough = [1:PCA_dim])
    
            weights = ICA.ica{m}.weights;
            iweights  = pinv(ICA.ica{m}.weights); % weights will be in output file from previous step
            finalics  = setdiff(1:PCA_dim,toremove{m}); % to
            TraMat{m} = iweights(:,finalics) * weights(finalics,:);
        end
        
        if isempty(cat(1,toremove{:}))
             aas_log(aap,false,'WARNING: No ICs to remove, but applying montage anyway');
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