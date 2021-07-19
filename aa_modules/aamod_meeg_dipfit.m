function [aap, resp] = aamod_meeg_dipfit(aap,task,subj,sess)

resp='';

switch task
    case 'report'
        for fn = cellstr(spm_select('FPList',aas_getsesspath(aap,subj,sess),'^diagnostic_.*jpg$'))'
            aap = aas_report_add(aap,subj,'<table><tr><td>');
            aap=aas_report_addimage(aap,subj,fn{1});
            aap = aas_report_add(aap,subj,'</td></tr></table>');
        end
    case 'doit'
        infname = cellstr(aas_getfiles_bystream(aap,'meeg_session',[subj sess],'meeg'));
        
        [~, EL] = aas_cache_get(aap,'eeglab');
        EL.load;
        
        indfnEEG = strcmp(spm_file(infname,'ext'),'set');
        EEG = pop_loadset('filepath',spm_file(infname{indfnEEG},'path'),'filename',spm_file(infname{indfnEEG},'filename'));
        
        volcondfile = aas_getsetting(aap,'volumeCondutionModel');
        if ~exist(volcondfile,'file'), volcondfile = fullfile(EL.dipfitPath,volcondfile); end
        if ~exist(volcondfile,'file'), aas_log(aap,true,sprintf('Volume condition model %s not found',aas_getsetting(aap,'volumeCondutionModel'))); end
        mrifile = aas_getfiles_bystream_multilevel(aap,'subject',subj,'structural');
        mri = ft_read_mri(mrifile,'dataformat','nifti_spm');
        mri.anatomy = mri.anatomy/max(mri.anatomy(:));
        
        trans = aas_getsetting(aap,'transformation');
        if ischar(trans) % target channel location
            if ~exist(trans,'file'), trans = fullfile(EL.dipfitPath,trans); end
            if ~exist(trans,'file'), aas_log(aap,true,sprintf('Channel location of the target %s not found',aas_getsetting(aap,'transformation'))); end
            [~,trans] = coregister(EEG.chanlocs, trans, 'warp', 'auto', 'manual', 'off');
        end
        
        % dipole settings
        EEG = pop_dipfit_settings( EEG,'chansel',1:EEG.nbchan, ...
            'hdmfile',volcondfile,'coordformat','MNI','mrifile',mrifile,'chanfile',aas_getfiles_bystream(aap,'channellayout'),...
            'coord_transform',trans ...            
            );
        % autofit: coarse + fine
        if isempty(EEG.icaweights) % no previuos ICA
            EEG = pop_multifit(EEG, ...
                'threshold',aas_getsetting(aap,'rejectionThreshold'),... % as recommended by Makoto ("otherwise you'll have trouble in making STUDY later")
                'dipplot','off',...
                'plotopt',{'normlen','on'});
        else % previuos ICA
            EEG = pop_multifit(EEG, 1:size(EEG.icaweights,1) ,...
                'threshold',aas_getsetting(aap,'rejectionThreshold'),... % as recommended by Makoto ("otherwise you'll have trouble in making STUDY later")
                'dipplot','off',...
                'plotopt',{'normlen','on'});
        end
        
        % symmetrically constrained bilateral dipoles
        if aas_getsetting(aap,'constrainSymmetrical'), EEG = fitTwoDipoles(EEG, 'LRR', 35); end
        
        % plot dipoles
        pop_dipplot( EEG, 1:numel(EEG.dipfit.model) ,'mri',mri,'normlen','on','view',[-1,-1,1],'gui','off');
        set(gcf,'PaperPositionMode','auto');
        print(gcf,'-noui',fullfile(aas_getsesspath(aap,subj,sess),['diagnostic_' mfilename '.jpg']),'-djpeg','-r150');
        close(gcf);
        
        % save
        outfname = spm_file(infname,'prefix','dipfit_');
        pop_saveset(EEG,'filepath',aas_getsesspath(aap,subj,sess),'filename',spm_file(outfname{indfnEEG},'basename'));

        EL.unload;
                
        %% Describe outputs
        aap = aas_desc_outputs(aap,subj,sess,'meeg',outfname);
    case 'checkrequirements'
        if ~aas_cache_get(aap,'eeglab'), aas_log(aap,true,'EEGLAB is not found'); end
end