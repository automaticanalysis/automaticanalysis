function [aap, resp] = aamod_meeg_converttoeeglab(aap,task,subj,sess)

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
        
        [junk, EL] = aas_cache_get(aap,'eeglab');
        EL.load;
        
        % read data
        EEG = [];
        for i = 1:numel(infname)
            try
                EEG = pop_fileio(infname{i});
                [res,o1,o2] = evalc('eeg_checkset(EEG)');
                if isempty(res), break;
                else, throw(res); end
            catch err
                E(i) = err;
            end
        end
        if isempty(EEG)
            for i = 1:numel(E)
                aas_log(aap,false,sprintf('ERROR: reading %s - %s',infname{i}, E(i).message));
            end
        end
        
        % channel layout
        EEG = pop_chanedit(EEG,'lookup',aas_getfiles_bystream(aap,'channellayout'));
        
        % remove channel
        if ~isempty(aas_getsetting(aap,'removechannel'))
            chns = strsplit(aas_getsetting(aap,'removechannel'),':');
            EEG = pop_select(EEG, 'nochannel', cellfun(@(x) find(strcmp({EEG.chanlocs.labels}, x)), chns));
        end
        
        % downsample
        if ~isempty(aas_getsetting(aap,'downsample'))
            sRate = aas_getsetting(aap,'downsample');
            if sRate ~= EEG.srate, EEG = pop_resample( EEG, aas_getsetting(aap,'downsample')); end
        end
                
        % diagnostics
        diagpath = fullfile(aas_getsesspath(aap,subj,sess),['diagnostic_' mfilename '_raw.jpg']);
        meeg_diagnostics_continuous(EEG,aas_getsetting(aap,'diagnostics'),'Raw',diagpath);
                
        fnameroot = sprintf('eeg_%s',aas_getsubjname(aap,subj));
        while ~isempty(spm_select('List',aas_getsesspath(aap,subj,sess),[fnameroot '.*']))
            fnameroot = ['aa' fnameroot];
        end
        
        pop_saveset(EEG,'filepath',aas_getsesspath(aap,subj,sess),'filename',fnameroot);
        outfname = spm_select('FPList',aas_getsesspath(aap,subj,sess),[fnameroot '.*']);
        
        EL.unload;
        
        %% Describe outputs
        aap = aas_desc_outputs(aap,subj,sess,'meeg',outfname);
    case 'checkrequirements'
        if ~aas_cache_get(aap,'eeglab'), aas_log(aap,true,'EEGLAB is not found'); end
end