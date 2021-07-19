function [aap, resp] = aamod_meeg_cleanartifacts(aap,task,subj,sess)

resp='';

switch task
    case 'report'
        for fn = cellstr(spm_select('FPList',aas_getsesspath(aap,subj,sess),'^diagnostic_.*jpg$'))'
            aap = aas_report_add(aap,subj,'<table><tr><td>');
            aap=aas_report_addimage(aap,subj,fn{1});
            aap = aas_report_add(aap,subj,'</td></tr></table>');
        end
        
        [~, EL] = aas_cache_get(aap,'eeglab');
        if subj == 1, EL.load;
        else, EL.reload; end
        outfname = cellstr(aas_getfiles_bystream(aap,'meeg_session',[subj sess],'meeg','output'));
        outfname = outfname{strcmp(spm_file(outfname,'ext'),'set')};
        EEG = pop_loadset('filepath',spm_file(outfname,'path'),'filename',spm_file(outfname,'filename'));
        EL.unload;
        
        if isfield(EEG.etc,'clean_channel_mask')
            keptData = mean(EEG.etc.clean_channel_mask)*mean(EEG.etc.clean_sample_mask);
        else
            keptData = mean(EEG.etc.clean_sample_mask);
        end
                
        aap = aas_report_add(aap,subj,sprintf('<p><b>Ratio of preserved data:</b> %3.2f%%</p>',keptData*100));
    case 'doit'
        infname = cellstr(aas_getfiles_bystream(aap,'meeg_session',[subj sess],'meeg'));
        
        [~, EL] = aas_cache_get(aap,'eeglab');
        EL.load;
        
        indfnEEG = strcmp(spm_file(infname,'ext'),'set');
        origEEG = pop_loadset('filepath',spm_file(infname{indfnEEG},'path'),'filename',spm_file(infname{indfnEEG},'filename'));
        params = {};
        criteria = aas_getsetting(aap,'criteria');
        for f = fieldnames(criteria)'
            if strcmp(f{1}, 'COMMENT'), continue; end
            if ~isempty(criteria.(f{1}))
                params{end+1} = f{1};
                params{end+1} = criteria.(f{1});
            end
        end
        params = [params {'MaxMem' aap.options.aaparallel.memory*1000*0.75}]; % 75% of the available memory
        
        EEG = clean_artifacts(origEEG,params{:});
        
        if isfield(EEG.etc,'clean_channel_mask')
            keptData = mean(EEG.etc.clean_channel_mask)*mean(EEG.etc.clean_sample_mask);
        else
            keptData = mean(EEG.etc.clean_sample_mask);
        end
        
        vis_artifacts(EEG,origEEG,'WindowLength',floor(EEG.pnts/EEG.srate));
        set(gcf,'position',[0,0,1920 600]);
        set(gcf,'PaperPositionMode','auto');
        set(gcf,'Name',sprintf('%2.3f%% of the data kept',keptData*100))
        print(gcf,'-noui',fullfile(aas_getsesspath(aap,subj,sess),['diagnostic_' mfilename '.jpg']),'-djpeg','-r150');
        close(gcf);
        
        if isfield(EEG.etc,'clean_channel_mask') && ~strcmp(aas_getsetting(aap,'interpolate'), 'off')
            EEG = pop_interp(EEG, origEEG.chanlocs, aas_getsetting(aap,'interpolate'));
        end
        
        outfname = spm_file(infname,'prefix','cleaned_');
        pop_saveset(EEG,'filepath',aas_getsesspath(aap,subj,sess),'filename',spm_file(outfname{indfnEEG},'basename'));
        
        EL.unload;
                
        %% Describe outputs
        aap = aas_desc_outputs(aap,subj,sess,'meeg',outfname);
    case 'checkrequirements'
        if ~aas_cache_get(aap,'eeglab'), aas_log(aap,true,'EEGLAB is not found'); end
end