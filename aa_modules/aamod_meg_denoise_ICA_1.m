function [aap, resp] = aamod_meg_denoise_ICA_1(aap,task,subj,sess)

resp='';

switch task
    case 'report'
        
    case 'doit'
        %% Initialise
        sessdir = aas_getsesspath(aap,subj,sess);
        infname = aas_getfiles_bystream(aap,'meg_session',[subj sess],'meg'); infname = basename(infname(1,:));
        outfname = fullfile(sessdir,[infname '_ICA']); % specifying output filestem
        
        %% Not sure if bit below is needed?
%         if aas_stream_has_contents(aap,[subj sess],'meg_ica')
%             ICA0 = load(aas_getfiles_bystream(aap,'meg_session',[subj sess],'meg_ica'));
%             ICA = struct(...
%                 'weights',ICA0.ica.weights, ...
%                 'compvars',ICA0.ica.compvars ...
%                 );
%         else
%             ICA=struct();    
%         end
               
        D = spm_eeg_load(fullfile(sessdir,infname));  % INPUTSTREAM from Convert

        modalities = textscan(aap.tasklist.currenttask.settings.modalities,'%s','delimiter',':'); modalities = modalities{1}';
         
        %% Could make these two parameters for AA
        PCA_dim = aap.tasklist.currenttask.settings.PCA_dim;           % Number of PCs for ICA
        Rseed = aap.tasklist.currenttask.settings.Rseed;               % to make reproducible
      
        %% RUN
        ica = {};
        for n = 1:numel(modalities)
            chans = indchantype(D,modalities{n});   
            if PCA_dim > length(chans)
                actual_PCA_dim = round(0.75*length(chans)); 
                aas_log(aap,false,'Warning: For modality %s, ICA PCA dimensions set to %d (because only %d sensors)',modalitiles{n},actual_PCA_dim,length(chans));
            else
                actual_PCA_dim = PCA_dim;
            end
            try
                if isempty(Rseed)
                    [weights,sphere,compvars,bias,signs,lrates,ICs] = runica(d,'pca',actual_PCA_dim,'extended',1,'maxsteps',800); % will give different answer each time run
                else
                    [weights,sphere,compvars,bias,signs,lrates,ICs] = rik_runica(D(chans,:),'pca',actual_PCA_dim,'extended',1,'maxsteps',800,'rseed',Rseed); % Just local copy where rand seed can be passed
                end
                ica{n}.weights = weights;
                ica{n}.chans = chans;
                ica{n}.compvars = compvars;
                ica{n}.modality = modalities{n};
            catch E
                aas_log(aap,true,['MEG:detect_ICA_artefacts:' E.message]);
            end
        end
        
        %% Outputs
        save(outfname,'ica');
        aap=aas_desc_outputs(aap,subj,sess,'meg_ica',[outfname '.mat']);
end