function [aap, resp] = aamod_meeg_rereference(aap,task,subj,sess)

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
        
        EEG = pop_loadset(infname{strcmp(spm_file(infname,'ext'),'set')});
        
        EEG.nbchan = EEG.nbchan+1;
        EEG.data(end+1,:) = zeros(1, EEG.pnts);
        EEG.chanlocs(EEG.nbchan).labels = 'initialReference';
        
        ref = aas_getsetting(aap,'reference');
        if ischar(ref)
            if strcmp(ref,'average'), ref = [];
            elseif any(ref==':'), ref = strsplit(ref,':');
            elseif ~any(strcmp({EEG.chanlocs(:,1).labels},ref))
                aas_log(aap,true,sprintf('Channel %s is not found',ref));
            end
        elseif isnumeric(ref)
            if max(ref) > size(EEG.chanlocs,1), aas_log(aap,true,sprintf('There is only %d channel',size(EEG.chanlocs,1))); end
        end
        
        EEG = pop_reref(EEG, ref);
        
        EEG = pop_select( EEG,'nochannel',{'initialReference'});
        
        % diagnostics
        diagpath = fullfile(aas_getsesspath(aap,subj,sess),['diagnostic_' mfilename '.jpg']);
        meeg_diagnostics_continuous(EEG,aas_getsetting(aap,'diagnostics'),'Re-referenced',diagpath);
        
        outfname = spm_file(infname,'prefix','reref_');
        pop_saveset(EEG,'filepath',aas_getsesspath(aap,subj,sess),'filename',spm_file(outfname{1},'basename'));
        
        EL.unload;
        
        %% Describe outputs
        aap = aas_desc_outputs(aap,subj,sess,'meeg',outfname);
    case 'checkrequirements'
        if ~aas_cache_get(aap,'eeglab'), aas_log(aap,true,'EEGLAB is not found'); end
        if isempty(aas_getsetting(aap,'reference')), aas_log(aap,true,'Reference MUST be specified.'); end
end