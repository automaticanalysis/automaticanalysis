function [aap, resp] = aamod_meg_denoise_ICA_1(aap,task,subj,sess)

resp='';

switch task
    case 'report'
        
    case 'doit'
        %% Initialise
        sessdir = aas_getsesspath(aap,subj,sess);
        infname = aas_getfiles_bystream(aap,'meg_session',[subj sess],'meg'); infname = basename(infname(1,:));
        outfname = fullfile(sessdir,[infname '_ICA']); % specifying output filestem
        if aas_stream_has_contents(aap,[subj sess],'meg_ica')
            ICA0 = load(aas_getfiles_bystream(aap,'meg_session',[subj sess],'meg_ica'));
            ICA = struct(...
                'weights',ICA0.ica.weights, ...
                'compvars',ICA0.ica.compvars ...
                );
        else
            ICA=struct();    
        end
        
        modalities = textscan(aap.tasklist.currenttask.settings.modalities,'%s','delimiter',':'); modalities = modalities{1}';
        arttopos=load(aas_getfiles_bystream(aap,'topography'));
        ref_chans = textscan(aap.tasklist.currenttask.settings.artifactdetection.ref_chans,'%s','delimiter',':'); ref_chans = ref_chans{1}';
        ref_type = textscan(aap.tasklist.currenttask.settings.artifactdetection.ref_type,'%s','delimiter',':'); ref_type = ref_type{1}';
        samp = aap.tasklist.currenttask.settings.sampling;
                
        D=spm_eeg_load(fullfile(sessdir,infname));  % INPUTSTREAM from Convert
        
        ICA.PCA_dim = 60;           % Number of PCs for ICA
        ICA.VarThr = 0;             % 100/PCA_dim;
        if aap.tasklist.currenttask.settings.artifactdetection.doPermutation >= 0 % -1 - automatic --> do not add field
            ICA.Nperm = aap.tasklist.currenttask.settings.artifactdetection.doPermutation;
        end
        ICA.FiltPars = [0.1 40];    % [0.1 40] %% [0.05 20];
        ICA.Rseed = 1;              % to make reproducible
        
        ICA.TemAbsPval = .05/ICA.PCA_dim; % 0.05;  % Too liberal if not permuted?
        ICA.TemRelZval = aap.tasklist.currenttask.settings.artifactdetection.TemRelZval;    % SETTINGS for detect artifact
        ICA.SpaAbsPval = .05/ICA.PCA_dim; % 0.05;  % Too liberal if not permuted?
        ICA.SpaRelZval = aap.tasklist.currenttask.settings.artifactdetection.SpaRelZval;    % SETTINGS for detect artifact % 3 is too strict (since high topo correlation anyway)
        ICA.thresholding = aap.tasklist.currenttask.settings.artifactdetection.thresholding;% SETTINGS for thresholding artifact
%        ICA.remove = aap.tasklist.currenttask.settings.artifactdetection.remove;            % SETTINGS for remove artifact
        
        ICA.refs.tem = {};  ICA.refs.spa = {};              % Reference signal for correlating with ICs
        if isempty(samp), samp = 1:size(D,2); end;
        for a = 1:length(ref_chans),
            ICA.refs.tem{a} = D(indchannel(D,ref_chans{a}),samp);
        end
        ICA.refs.type = ref_type;
        
        %% RUN
        chans = {}; 
        for n = 1:numel(modalities)
            m = find(strcmp(arttopos.modalities,modalities{n}));
            for t = 1:numel(ICA.refs.type)
                tmp = getfield(arttopos,ICA.refs.type{t}); 
                ICA.refs.spa{t} = tmp{m}(:)'; % detect_ICA_artefacts below assumes 1xChan vector
            end
            %if ~aap.tasklist.currenttask.settings.artifactdetection.useECGtopo, ICA.refs.spa = []; end  % ECG topography may be too variable
            chans{n} = find(strcmp(D.chantype,modalities{n}));  
            ICA.d  = D(chans{n},:);
            ICA.FiltPars = [ICA.FiltPars D.fsample];
            try
                ica{n} = detect_ICA_artefacts(ICA);
            catch E
                aas_log(aap,true,['MEG:detect_ICA_artefacts:' E.message]);
            end
        end
        
        %% Outputs
        save(outfname,'chans','ica');
        aap=aas_desc_outputs(aap,subj,sess,'meg_ica',[outfname '.mat']);
end