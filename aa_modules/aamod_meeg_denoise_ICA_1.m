function [aap, resp] = aamod_meeg_denoise_ICA_1(aap,task,subj,sess)

resp='';

switch task
    case 'report'
        
    case 'doit'
        [s,EL] = aas_cache_get(aap,'eeglab');
        if ~s, aas_log(aap,true,'EEGLAB not found'); end
        EL.load;
        
        %% Initialise
        [junk, SPMtool] = aas_cache_get(aap,'spm');
        SPMtool.doToolbox('fieldtrip','load');
        SPMtool.addCollection('meeg');
        [s,EL] = aas_cache_get(aap,'eeglab');
        if ~s, aas_log(aap,true,'EEGLAB not found'); end
        EL.load;
        
        sessdir = aas_getsesspath(aap,subj,sess);
        infname = aas_getfiles_bystream(aap,'meeg_session',[subj sess],'meeg'); infname = basename(infname(1,:));
        outfname = fullfile(sessdir,[infname '_ICA']); % specifying output filestem
        
        %% Not sure if bit below is needed?
%         if aas_stream_has_contents(aap,[subj sess],'meg_ica')
%             ICA0 = load(aas_getfiles_bystream(aap,'meeg_session',[subj sess],'meg_ica'));
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
%                 if ~isempty(Rseed) && exist('rik_runica','file')              
%                     [weights,sphere,compvars,bias,signs,lrates,ICs] = rik_runica(D(chans,:),'pca',actual_PCA_dim,'extended',1,'maxsteps',800,'rseed',Rseed); % Just local copy where rand seed can be passed
%                 else
                    if ~isempty(Rseed), aas_log(aap,false,'WARNING: Random seed requested but no facility with standard runica?'); end
                    [weights,sphere,compvars,bias,signs,lrates,ICs] = runica(D(chans,:),'pca',actual_PCA_dim,'extended',1,'maxsteps',800); % will give different answer each time run
%                 end
                ica{n}.weights = weights;
                ica{n}.chans = chans;
                ica{n}.compvars = compvars;
                ica{n}.modality = modalities{n};
            catch E
                aas_log(aap,true,['MEG:detect_ICA_artefacts:' E.message]);
            end
        end
        
        EL.unload;
        SPMtool.rmCollection('meeg');
        SPMtool.doToolbox('fieldtrip','unload')

        %% Outputs
        save(outfname,'ica');
        aap=aas_desc_outputs(aap,subj,sess,'meeg_ica',[outfname '.mat']);
    case 'checkrequirements'
        if ~aas_cache_get(aap,'spm'), aas_log(aap,true,'SPM is not found'); end
end