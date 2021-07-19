function [aap, resp] = aamod_meeg_ica(aap,task,subj,sess)

resp='';

switch task
    case 'report'
        
    case 'doit'
        infname = cellstr(aas_getfiles_bystream(aap,'meeg_session',[subj sess],'meeg'));
        
        [~, EL] = aas_cache_get(aap,'eeglab');
        EL.load;
        
        indfnEEG = strcmp(spm_file(infname,'ext'),'set');
        EEG = pop_loadset('filepath',spm_file(infname{indfnEEG},'path'),'filename',spm_file(infname{indfnEEG},'filename'));
        
        % common parameters
        if isempty(aas_getsetting(aap,'PCA')), PCA = 0;
        elseif isnumeric(aas_getsetting(aap,'PCA')), PCA = aas_getsetting(aap,'PCA');
        elseif ischar(aas_getsetting(aap,'PCA')) && strcmp(aas_getsetting(aap,'PCA'),'rank')
            if isfield(EEG.etc, 'clean_channel_mask')
                PCA = min([rank(double(EEG.data')) sum(EEG.etc.clean_channel_mask)]);
            else
                PCA = rank(double(EEG.data'));
            end
        else, aas_log(aap,true,'Unknown PCA option');
        end
        
        iter = 0;
        if ~isempty(aas_getsetting(aap,'iterations')), iter = aas_getsetting(aap,'iterations'); end
        
        params = {};
        switch aas_getsetting(aap,'method')
            case 'AMICA'
                if PCA, params = {'pcakeep',PCA}; end
                if iter, params = [params {'max_iter' iter}]; end
                if ~isempty(aas_getsetting(aap,'options.AMICA.num_models')), params = [params {'num_models' aas_getsetting(aap,'options.AMICA.num_models')}]; end
                if ~isempty(aas_getsetting(aap,'options.AMICA.numrej'))
                    params = [params {'do_reject' 1 'numrej' aas_getsetting(aap,'options.AMICA.numrej')}]; 
                    if ~isempty(aas_getsetting(aap,'options.AMICA.rejint')), params = [params {'rejint' aas_getsetting(aap,'options.AMICA.rejint')}]; end
                    if ~isempty(aas_getsetting(aap,'options.AMICA.rejsig')), params = [params {'rejsig' aas_getsetting(aap,'options.AMICA.rejsig')}]; end
                end
                
                tmpdir = tempname;
                mkdir(tmpdir);
                currdir = pwd;
                cd(tmpdir);
                runamica15(EEG.data,'outdir',tmpdir,params{:});
                cd(currdir);
                EEG.etc.amica  = loadmodout15(tmpdir);
                EEG.etc.amica.S = EEG.etc.amica.S(1:EEG.etc.amica.num_pcs, :);
                EEG.icaweights = EEG.etc.amica.W;
                EEG.icasphere  = EEG.etc.amica.S;
                rmdir(tmpdir,'s');
                
            case 'runica'
                if PCA, params = {'pca',PCA}; end
                if ~isempty(aas_getsetting(aap,'options.runica.extended')), params = [params {'extended' aas_getsetting(aap,'options.runica.extended')}]; end
                if iter, params = [params {'maxsteps' iter}]; end
                
                EEG = pop_runica(EEG,'icatype','runica',params{:});
        end

        % check output
        [res,EEG] = evalc('eeg_checkset(EEG,''ica'')');
        if ~isempty(res) && aas_getsetting(aap,'errorOnFailedCheck'), aas_log(aap,true,res); end
        
        % save
        outfname = spm_file(infname,'prefix','ica_');
        pop_saveset(EEG,'filepath',aas_getsesspath(aap,subj,sess),'filename',spm_file(outfname{indfnEEG},'basename'));

        EL.unload;
                
        %% Describe outputs
        aap = aas_desc_outputs(aap,subj,sess,'meeg',outfname);
    case 'checkrequirements'
        if ~aas_cache_get(aap,'eeglab'), aas_log(aap,true,'EEGLAB is not found'); end
end