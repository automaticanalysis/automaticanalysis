function [aap, resp] = aamod_meg_denoise_ICA_1(aap,task,subj,sess)

resp='';

switch task
    case 'report'
        
    case 'doit'
        %% Initialise
        sessdir = aas_getsesspath(aap,subj,sess);
        infname = aas_getfiles_bystream(aap,'meg_session',[subj sess],'meg'); infname = basename(infname(1,:));
        outfname = fullfile(sessdir,[infname '_ICA']); % specifying output filestem
        
        modalities = textscan(aap.tasklist.currenttask.settings.modalities,'%s','delimiter',':'); modalities = modalities{1}';
        arttopos=load(aas_getfiles_bystream(aap,'topography'));
        ref_chans = textscan(aap.tasklist.currenttask.settings.artifactdetection.ref_chans,'%s','delimiter',':'); ref_chans = ref_chans{1}';
        samp = aap.tasklist.currenttask.settings.sampling;
                
        D=spm_eeg_load(fullfile(sessdir,infname));  % INPUTSTREAM from Convert
        
        ICA=struct();
        ICA.PCA_dim = 60;           % Number of PCs for ICA
        ICA.VarThr = 0;             % 100/PCA_dim;
        ICA.Nperm = 0;              % round((4*sqrt(CorPval)/CorPval)^2); % If want accurate p-values
        ICA.FiltPars = [0.1 40];    %[1 40] %% [0.05 20];
        ICA.Rseed = 1;              % to make reproducible
        
        ICA.TemAbsPval = .05/ICA.PCA_dim; % 0.05;  % Too liberal if not permuted?
        ICA.TemRelZval = aap.tasklist.currenttask.settings.artifactdetection.TemRelZval;         % SETTINGS for detect artifact
        ICA.SpaAbsPval = .05/ICA.PCA_dim; % 0.05;  % Too liberal if not permuted?
        ICA.SpaRelZval = aap.tasklist.currenttask.settings.artifactdetection.SpaRelZval;         % SETTINGS for detect artifact % 3 is too strict (since high topo correlation anyway)
        
        ICA.refs.tem = {};                % Reference signal for correlating with ICs
        if isempty(samp), samp = 1:size(D,2); end;
        for a = 1:length(ref_chans),
            ICA.refs.tem{a} = D(find(strcmp(D.chanlabels,ref_chans{a})),samp);
        end
        
        %% RUN
        chans = {}; remove = {}; weights = {}; temcor = {}; spacor = {}; TraMat = {};
        for m = 1:numel(modalities)
            ICA.refs.spa = {arttopos.HEOG{m}', arttopos.VEOG{m}', arttopos.ECG{m}'};  % Assumes modalities ordered same way!!!
            if m == 2, ICA.refs.spa{3} = []; end  % ECG topography too variable, so just select on timecourse?
            chans{m} = find(strcmp(D.chantype,modalities{m}));  % RH: Bad coding: modalities{2} may not be Grads - see other comment about only converting grads
            ICA.d  = D(chans{m},:);
            ICA.FiltPars = [ICA.FiltPars D.fsample];
            [remove{m},weights{m},TraMat{m},temcor{m},spacor{m}] = detect_ICA_artefacts(ICA);
        end
        
        %% Outputs
        save(outfname,'remove','TraMat','temcor','spacor','weights','chans');
        aap=aas_desc_outputs(aap,subj,sess,'meg_ica',[outfname '.mat']);
end