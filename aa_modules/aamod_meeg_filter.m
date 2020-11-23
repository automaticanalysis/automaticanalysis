function [aap, resp] = aamod_meeg_filter(aap,task,subj,sess)

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
        
        %% Load toolboxes
        [junk, FT] = aas_cache_get(aap,'fieldtrip');
        [junk, EL] = aas_cache_get(aap,'eeglab');

        FT.load;
        
        %% Run
        cfg = [];
        cfg.dataset = infname{strcmp(spm_file(infname,'ext'),'set')};
        
        for filt = {'lp','hp','bp','bs','dft'}
           conf = aas_getsetting(aap,[filt{1} 'freq']);
           if ~isempty(conf)
               cfg.([filt{1} 'filter']) = 'yes';
               cfg.([filt{1} 'freq']) = conf;
           end
        end
        if isfield(cfg,'dftfilter')
            cfg.dftreplace = aas_getsetting(aap,'dftreplace');
            cfg.dftbandwidth = aas_getsetting(aap,'dftbandwidth');
            if numel(cfg.dftbandwidth) < numel(cfg.dftfreq), cfg.dftbandwidth(end+1:numel(cfg.dftfreq)) = cfg.dftbandwidth(end); end
            cfg.dftneighbourwidth = aas_getsetting(aap,'dftneighbourwidth');
            if numel(cfg.dftneighbourwidth) < numel(cfg.dftfreq), cfg.dftneighbourwidth(end+1:numel(cfg.dftfreq)) = cfg.dftneighbourwidth(end); end
        end        
        if aas_getsetting(aap,'medianfilter'), cfg.medianfilter = 'yes'; end

        dat = ft_preprocessing(cfg);
        events = ft_read_event(cfg.dataset);
        
        FT.unload;
        rmpath(fullfile(FT.toolPath,'external','eeglab'));
        
        %% Load into EEGLAB
        EL.load;
        
        EEG = pop_fileio(dat.hdr, dat.trial{1}, events);
                
        % diagnostics
        diagpath = fullfile(aas_getsesspath(aap,subj,sess),['diagnostic_' mfilename '_filtered.jpg']);
        meeg_diagnostics_continuous(EEG,aas_getsetting(aap,'diagnostics'),'Filtered',diagpath);
        
        outfname = spm_file(infname,'prefix','filtered_');
        pop_saveset(EEG,'filepath',aas_getsesspath(aap,subj,sess),'filename',spm_file(outfname{1},'basename'));
        
        EL.unload;
        
        %% Describe outputs
        aap = aas_desc_outputs(aap,subj,sess,'meeg',outfname);
    case 'checkrequirements'
        if ~aas_cache_get(aap,'eeglab'), aas_log(aap,true,'EEGLAB is not found'); end
        if ~aas_cache_get(aap,'fieldtrip'), aas_log(aap,true,'FieldTrip is not found'); end
end