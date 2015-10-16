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
       
        samp = aap.tasklist.currenttask.settings.sampling;
        if isempty(samp),
            samp = [1:D.nsamples];
        end
        if rem(length(samp),2)==0  % Must be odd number of samples any fft below
            samp(1) = [];
        end
        
        arttopos = load(aas_getfiles_bystream(aap,'topography'));
        ref_chans = textscan(aap.tasklist.currenttask.settings.ref_chans,'%s','delimiter',':'); ref_chans = ref_chans{1}';
        ref_type = textscan(aap.tasklist.currenttask.settings.ref_type,'%s','delimiter',':'); ref_type = ref_type{1}';

        refs.tem = {};  % Temporal reference signal for correlating with ICs
        for r = 1:length(ref_chans),
            refs.tem{r} = D(indchannel(D,ref_chans{r}),samp);
        end

        thresholding = aap.tasklist.currenttask.settings.thresholding;    % SETTINGS for thresholding artifact
        if isempty(intersect(thresholding,[1:4]))
            aas_log(aap,true,'ERROR: thresholding option must be one or more of [1:4]');
        end
        
        % Not sure why just can't set Nperm = aap.tasklist.currenttask.settings.doPermutation... ?
        if aap.tasklist.currenttask.settings.doPermutation >= 0 % -1 - automatic --> do not add field
            Nperm = aap.tasklist.currenttask.settings.doPermutation;
        else
            Nperm = 0;
        end
        FiltPars = [0.1 40 D.fsample];    % [0.1 40] %% [0.05 20];
        
        PCA_dim = size(ICA.ica{1}.weights,1)
        
        TemAbsPval = aap.tasklist.currenttask.settings.TemAbsPval;    % SETTINGS for detect artifact
        TemRelZval = aap.tasklist.currenttask.settings.TemRelZval;    % SETTINGS for detect artifact
        SpaAbsPval = aap.tasklist.currenttask.settings.SpaAbsPval;    % SETTINGS for detect artifact
        SpaRelZval = aap.tasklist.currenttask.settings.SpaRelZval;    % SETTINGS for detect artifact
        VarThr = aap.tasklist.currenttask.settings.VarThr;            % SETTINGS for detect artifact
         
        %% Main loop
        toremove = cell(1,numel(ICA.ica)); TraMat = cell(1,numel(ICA.ica));       
        for m = 1:numel(ICA.ica)
            
            weights   = ICA.ica{m}.weights;
            iweights  = pinv(weights); % weights will be in output file from previous step
            
            ICs = weights * D(ICA.ica{m}.chans,samp);
            
            n = find(strcmp(arttopos.modalities,ICA.ica{m}.modality));
            refs.spa = {};
            for r = 1:numel(ref_type)
                tmp = getfield(arttopos,ref_type{r});
                refs.spa{r} = tmp{n}(:)'; % detect_ICA_artefacts below assumes 1xChan vector
            end
            
            %% Filtering (if any) (and transposition for speed)
            if length(FiltPars) == 3
                fprintf('Bandpass filtering from %d to %d Hz (warning - can fail)\n',FiltPars(1), FiltPars(2));
                for r=1:length(refs.tem)
                    refs.tem{r} = ft_preproc_bandpassfilter(refs.tem{r}, FiltPars(3), FiltPars(1:2),  [], 'but','twopass','reduce');
                end
                ICs  = ft_preproc_bandpassfilter(ICs,  FiltPars(3), FiltPars(1:2),  [], 'but','twopass','reduce')';
            elseif length(FiltPars) == 2
                fprintf('Lowpass filtering to %d Hz\n',FiltPars(1));
                for r=1:length(refs.tem)
                    refs.tem{r} = ft_preproc_lowpassfilter(refs.tem{r}, FiltPars(2), FiltPars(1),  5, 'but','twopass','reduce');
                end
                ICs  = ft_preproc_lowpassfilter(ICs,  FiltPars(2), FiltPars(1),  5, 'but','twopass','reduce')';
            else
                ICs  = ICs'; % Faster for FFT  if pre-transpose once (check tic;toc)
            end
            % figure; for i=1:Nrefs; plot(refs.tem{r}'); title(i); pause; end
            % figure; for i=1:PCA_dim; plot(ICs(:,i)); title(i); pause; end
            % figure; for i=1:PCA_dim; plot(ICs(30000:40000,i)); title(i); pause; end
            
            tempcor = zeros(length(refs.tem),PCA_dim); tempval = zeros(length(refs.tem),PCA_dim);           
            temprem  = cell(4,length(refs.tem));          
            for r = 1:length(refs.tem)
                if ~isempty(refs.tem{r})    % Check temporal correlation with any reference channels
                    
                    for k = 1:PCA_dim
                        [tempcor(r,k),tempval(r,k)] = corr(refs.tem{r}',ICs(:,k));
                    end
                    
                    [~,temprem{1,r}] = max(abs(tempcor(r,:)));
                    
                    temprem{2,r} = find(tempval(r,:) < TemAbsPval);
                    
                    temprem{3,r} = find(abs(zscore(tempcor(r,:))) > TemRelZval);
                    
                    if Nperm > 0
                        permcor = zeros(1,PCA_dim);
                        maxcor  = zeros(Nperm,1);
                        
                        ff = fft(refs.tem{r}',Nsamp);
                        mf = abs(ff);
                        wf = angle(ff);
                        hf = floor((length(ff)-1)/2);
                        rf = mf;
                        
                        for l = 1:Nperm % could parfor...
                            rf(2:hf+1)=mf(2:hf+1).*exp((0+1i)*wf(randperm(hf)));    % randomising phases (preserve mean, ie rf(1))
                            rf((hf+2):length(ff))=conj(rf((hf+1):-1:2));            % taking complex conjugate
                            btdata = ifft(rf,Nsamp);                                % Inverse Fourier transform
                            
                            for k = 1:PCA_dim
                                permcor(k) = corr(btdata,ICs(:,k));
                            end
                            maxcor(l) = max(abs(permcor));
                            fprintf('.');
                        end
                        fprintf('\n')
                        %         figure,hist(maxcor)
                        
                        temprem{4,r} = find(abs(tempcor(r,:)) > prctile(maxcor,100*(1-PermPval)));
                    end
                else
                    temprem{1,r} = []; temprem{2,r} = []; temprem{3,r} = []; temprem{4,r} = [];
                end
            end
            
            spatcor = zeros(length(refs.spa),PCA_dim); spatval = zeros(length(refs.spa),PCA_dim);
            spatrem  = cell(3,length(refs.spa));
            for r = 1:length(refs.spa)
                if ~isempty(refs.spa{r})
                    %% Check spatial correlation with any reference channels
                    for k = 1:PCA_dim
                        [spatcor(r,k),spatval(r,k)] = corr(refs.spa{r}',iweights(:,k));
                    end
                    
                    [~,spatrem{1,r}] = max(abs(spatcor(r,:)));
                    
                    spatrem{2,r} = find(spatval(r,:) < SpaAbsPval);
                    
                    spatrem{3,r} = find(abs(zscore(spatcor(r,:))) > SpaRelZval);
                else
                    spatrem{1,r} = []; spatrem{2,r} = []; spatrem{3,r} = []; spatrem{4,r} = [];
                end
            end

            if VarThr > 0 % Variance Thresholding
                compvars  = ICA.ica{m}.compvars;
                varexpl   = 100*compvars/sum(compvars);
                varenough = find(varexpl > VarThr);
            else
                varenough = [1:PCA_dim];
            end
                      
            %% Define toremove                  
            toremove{m} = [1:PCA_dim];
            switch aap.tasklist.currenttask.settings.toremove
                case 'temp' % define artefacts based on temporal correlation thresholds only
                    for t = 1:length(thresholding)
                        toremove{m} = intersect(toremove{m},unique(cat(2,temprem{thresholding(t),:})));
                    end
                case 'spat' % define artefacts based on spatial correlation thresholds only
                    for t = 1:length(thresholding)
                        toremove{m} = intersect(toremove{m},unique(cat(2,spatrem{thresholding(t),:})));
                    end
                case 'either' % define artefacts based on either spatial OR temporal correlation thresholds
                    tempremove = [1:PCA_dim];
                    for t = 1:length(thresholding)
                        tempremove = intersect(tempremove,unique(cat(2,temprem{thresholding(t),:})));
                    end
                    spatremove = [1:PCA_dim];
                    for t = 1:length(thresholding)
                        spatremove = intersect(spatremove,unique(cat(2,spatrem{thresholding(t),:})));
                    end
                    toremove{m} = unique([tempremove spatremove]);
                case 'both' % define artefacts that meet both temporal AND spatial correlation thresholds (DEFAULT)
                    % !!!NOTE: This assumes temporal and spatial references were paired!!!
                    if size(temprem,2) ~= size(spatrem,2)
                        aas_log(aap,true,'ERROR: Must have same number of temporal and spatial references to use "both" option');
                    else
                        toremove{m} = [];
                        for r = 1:size(temprem,2)
                            bothremove = [1:PCA_dim];
                            for t = 1:length(thresholding)  % "If" bit below to handle cases where no reference for some modalities (eg EEG topo), even if exist for other modalities
                                if isempty(refs.tem{r}) & ~isempty(refs.spa{r})
                                    both = spatrem{thresholding(t),r};
                                elseif ~isempty(refs.tem{r}) & isempty(refs.spa{r})
                                    both = temprem{thresholding(t),r};
                                else
                                    both = intersect(temprem{thresholding(t),r},spatrem{thresholding(t),r});
                                end
                                bothremove = intersect(bothremove,both);
                            end
                            toremove{m} = [toremove{m} bothremove];
                        end
                    end
            end
            
            toremove{m} = intersect(toremove{m},varenough); % Additional criterion of sufficient variance (normally ignored, ie varenough = [1:PCA_dim])
            
            finalics  = setdiff(1:PCA_dim,toremove{m}); % to
            TraMat{m} = iweights(:,finalics) * weights(finalics,:);
        end
               
        if isempty(cat(2,toremove{:}))
            aas_log(aap,false,'WARNING: No ICs to remove - but applying identity montage anyway');
        end
 
        %% RUN
        S = []; S.D = D; achans = [];
        for m = 1:numel(ICA.ica)
            achans = [achans ICA.ica{m}.chans];
        end
        for c = 1:length(achans)
            S.montage.labelorg{c} = D.chanlabels{achans(c)};
        end
        S.montage.labelnew = S.montage.labelorg;
        S.montage.tra      = blkdiag(TraMat{:});
        S.keepothers       = 1;
        S.prefix           = 'M';
        D = spm_eeg_montage(S);
        
        %% Outputs
        sessdir = aas_getsesspath(aap,subj,sess);
        infname = aas_getfiles_bystream(aap,'meg_session',[subj sess],'meg'); infname = basename(infname(1,:));
        outfname = fullfile(sessdir,[infname '_ICA']); % specifying output filestem
        save(outfname,'toremove','TraMat','tempval','tempcor','spatval','spatcor','varenough');
        
        outfname = [S.prefix basename(infname_meg)];
        aap=aas_desc_outputs(aap,subj,sess,'meg',char([outfname '.dat'],[outfname '.mat']));
end