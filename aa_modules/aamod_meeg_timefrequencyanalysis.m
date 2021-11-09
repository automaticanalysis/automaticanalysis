function [aap, resp] = aamod_meeg_timefrequencyanalysis(aap,task,subj)

resp='';

switch task
    case 'report'
        bands = aas_getsetting(aap,'diagnostics.snapshotfwoi'); bands = ['multiplot'; mat2cell(bands,ones(1,size(bands,1)))];
        models = strrep(spm_file(cellstr(aas_getfiles_bystream(aap,'subject',subj,'timefreq')),'basename'),'timefreq_','')';
        aap = aas_report_add(aap,subj,'<table id="data"><tr>');
        aap = aas_report_add(aap,subj,'<th>Bands</th>');
        for m = models, aap = aas_report_add(aap,subj,['<th>Model: ' m{1} '</th>']); end
        aap = aas_report_add(aap,subj,'</tr>');
        
        for b = bands'
            aap = aas_report_add(aap,subj,'<tr>');
            if ischar(b{1})
                aap = aas_report_add(aap,subj,['<td>' b{1} '</td>']);
            else
                aap = aas_report_add(aap,subj,sprintf('<td>%1.2f-%1.2f</td>',b{1}));
            end
            for m = models
                aap = aas_report_add(aap,subj,'<td>');
                if ischar(b{1})
                    fnimg = fullfile(aas_getsubjpath(aap,subj),['diagnostic_' mfilename  '_' m{1} '_' b{1} '.jpg']);
                else
                    fnimg = fullfile(aas_getsubjpath(aap,subj),['diagnostic_' mfilename  '_' m{1} '_topoplot_freq-' sprintf('%1.2f-%1.2f',b{1}) '.jpg']);
                end
                if exist(fnimg,'file'), aap=aas_report_addimage(aap,subj,fnimg); end
                aap = aas_report_add(aap,subj,'</td>');
            end
            aap = aas_report_add(aap,subj,'</tr>');
        end
        
        aap = aas_report_add(aap,subj,'</table>');
    case 'doit'
        [~, FT] = aas_cache_get(aap,'fieldtrip');
        FT.load;
        FT.addExternal('spm12');
        
        tfa = aas_getsetting(aap,'timefrequencyanalysis');
        tfacfg = keepfields(tfa,{'method','taper','foi'});
        tfacfg.precision = 'single';
        tfacfg.pad         = 'nextpow2';
        tfacfg.output      = 'powandcsd';
        if ~isempty(tfa.twoicps), tfacfg.t_ftimwin = tfa.twoicps./tfacfg.foi; end
        if ~isempty(tfa.spectralsmoothing), tfacfg.tapsmofrq = tfa.spectralsmoothing*tfacfg.foi; end
        if ~isempty(tfa.toi) && isnumeric(tfa.toi), tfacfg.toi = tfa.toi/1000; end        
        tfacfg.keeptrials  = 'no';
        
        % baseline correction
        baswin = aas_getsetting(aap,'baselinewindow');
        if ~isempty(baswin)
            if isnumeric(baswin), baswin = baswin/1000; end % convert to seconds
        end
        bccfg = [];
        bccfg.baseline     = baswin;
        bccfg.baselinetype = 'relative';
        bccfg.parameter    = {'powspctrm' 'crsspctrm'};
        
        diag = aas_getsetting(aap,'diagnostics');
        diag.parameter = 'powspctrm';
        
        models = aas_getsetting(aap,'trialmodel');
        subjmatches=strcmp(aap.acq_details.subjects(subj).subjname,{models.subject});
        if ~any(subjmatches), aas_log(aap,true,['No trialmodel specification found for ' aas_getsubjdesc(aap,subj)]); end
        
        combinecfg = [];
        combinecfg.parameter = {'powspctrm' 'crsspctrm'};
        combinecfg.normalise = 'no';
        
        for m = models(subjmatches).model
            
            conSessions = cell(1,numel(m.session.names));
            
            % process sessions
            includedsessionnumbers = cellfun(@(x) find(strcmp({aap.acq_details.meeg_sessions.name},x)),m.session.names);
            for sess = 1:numel(includedsessionnumbers)
                sessnum = (includedsessionnumbers(sess));
                
                clear data
                meegfn = cellstr(aas_getfiles_bystream(aap,'meeg_session',[subj sessnum],'meeg'));
                switch spm_file(meegfn{1},'ext')
                    case 'mat'
                        filetype = 'fieldtrip';
                    case 'set' 
                        filetype = 'eeglab';
                        meegfn = meegfn(strcmp(spm_file(meegfn,'ext'),'set'));
                    otherwise
                        aas_log(aap,true,'Unsupported file format')
                end
                for seg = 1:numel(meegfn)
                    switch filetype
                        case 'fieldtrip'
                            dat = load(meegfn{seg});
                            data(seg) = ft_struct2single(dat.data);
                        case 'eeglab'
                            FT.unload;
                            if seg == 1
                                [~, EL] = aas_cache_get(aap,'eeglab');
                                EL.load;
                            else
                                EL.reload;
                            end
                            EEG = pop_loadset('filepath',spm_file(meegfn{seg},'path'),'filename',spm_file(meegfn{seg},'filename'));
                            if isempty(EEG.epoch)
                                aas_log(aap,false,sprintf('WARNING: segment # %d has no trial --> skipped',seg));
                                continue; 
                            end
                            EL.unload;
                            FT.reload;                            
                            data(seg) = ft_struct2single(eeglab2fieldtripER(EEG,'reorient',1));                            
                    end
                    
                    % select data
                    if isfield(data(seg),'ureventinfo')
                        toremove = [];
                        if ~isempty(aas_getsetting(aap,'ignorebefore'))
                            lim = aas_getsetting(aap,'ignorebefore');
                            field = 'eventnum_all';
                            if lim < 0
                                lim = -lim;
                                field = 'eventnum';
                            end
                            toremove = data(seg).ureventinfo.(field) < lim;
                        end
                        if ~isempty(aas_getsetting(aap,'ignoreafter'))
                            lim = aas_getsetting(aap,'ignoreafter');
                            field = 'eventnum_all';
                            if lim < 0
                                lim = -lim;
                                field = 'eventnum';
                            end
                            toremove = data(seg).ureventinfo.(field) > lim;
                        end
                        
                        if any(toremove)
                            cfg = [];
                            cfg.trials = ~toremove;
                            data(seg) = ft_selectdata(cfg,data(seg));
                        end
                    else
                        aas_log(aap,false,'WARNING: original eventinfo (ureventinfo) is not available -> ignorebefore and ignoreafter will be ignored')
                    end
                end
                data(cellfun(@isempty, {data.trial})) = []; % remove skipped segments
                
                % process events
                kvs = regexp(spm_file(meegfn{1},'basename'),'[A-Z]+-[0-9]+','match');
                events = cellfun(@(x) strsplit(x,'-'), kvs,'UniformOutput',false);
                conEvents = cell(1,numel(m.event.names));
                
                for e = 1:numel(m.event.names)
                    eventLabel = m.event.names{e};
                    trialinfo = str2double(events{cellfun(@(x) strcmp(x{1},eventLabel), events)}{2});
                    aas_log(aap,false,sprintf('INFO: processing event %s with trialinfo %d',eventLabel,trialinfo));
                    
                    % main 
                    tf = {}; weights = [];
                    for i = 1:numel(data)
                        cfg = tfacfg;
                        cfg.trials = find(data(i).trialinfo==trialinfo);
                        if isempty(cfg.trials), continue; end
                        if isfield(cfg,'toi') && ischar(cfg.toi) && strcmp(cfg.toi,'all') % whole trial
                            cfg.toi = (data(1).time{1}(1)+data(1).time{1}(end))/2; % centre
                            cfg.t_ftimwin = (data(i).time{1}(end)-data(i).time{1}(1))*ones(1,numel(cfg.foi));
                        end
                        tf{end+1} = ft_freqanalysis(cfg, data(i));
                        % baseline correction
                        if ~isempty(baswin)
                            lbc = tf{end}.labelcmb;
                            tf{end} = ft_freqbaseline(bccfg,tf{end}); 
                            tf{end}.labelcmb = lbc;
                        end
                        if aas_getsetting(aap,'weightedaveraging')
                            weights(end+1) = numel(cfg.trials);
                        else
                            weights(end+1) = 1;
                        end
                    end
                    if isempty(tf), continue; end
                    
                    cfg = combinecfg; 
                    cfg.normalise = 'yes';
                    cfg.weights = weights;
                    timefreqMain = ft_combine(cfg,tf{:});
                    
                    diagFn = fullfile(aas_getsesspath(aap,subj,sess),['diagnostic_' mfilename '_' eventLabel]);
                    if ~(ischar(m.samplevector) && strcmp(m.samplevector,'cont')) &&... % not for continuous 
                        ~exist([diagFn '_multiplot.jpg'],'file')
                        meeg_diagnostics_TFR(timefreqMain,diag,eventLabel,diagFn);
                    end
                    
                    % trialmodel
                    timefreqModel = timefreqMain;
                    if ischar(m.samplevector)
                        switch m.samplevector
                            case 'avg' % average - it is done
                            % do nothing
                            case 'cont' % continuously sampled data
                                clear tf
                                for i = 1:numel(data)
                                    cfg = keepfields(tfacfg,{'precision','pad','output','foi'});
                                    cfg.method = 'mtmfft';
                                    cfg.taper = 'hanning';
                                    cfg.trials = find(data(i).trialinfo==trialinfo);
                                    cfg.keeptrials = 'yes';
                                    tmptfr = ft_freqanalysis(cfg, data(i));
                                    
                                    % convert trials to time
                                    tf{i}           = tmptfr;
                                    tf{i}.powspctrm = permute(tmptfr.powspctrm, [2, 3, 1]);
                                    tf{i}.crsspctrm = permute(tmptfr.crsspctrm, [2, 3, 1]);
                                    tf{i}.dimord    = 'chan_freq_time'; % it used to be 'rpt_chan_freq'
                                    if isfield(data(i),'ureventinfo')
                                        tf{i}.time      = data(i).ureventinfo.latency(cfg.trials);
                                    else
                                        aas_log(aap,false,'WARNING: original eventinfo (ureventinfo) is not available -> sampleinfo will be used')
                                        tf{i}.time = data(i).sampleinfo(cfg.trials,1)/data(i).fsample;
                                    end
                                end
                                cfg = combinecfg;
                                cfg.normalise = 'yes';
                                cfg.weights = weights;
                                timefreqModel = ft_combine(cfg,tf{:});
                            case 'segmentavg'
                                timefreqModel.dimord    = 'chan_freq_time';
                                if isfield(data(1),'ureventinfo')
                                    timefreqModel.time = arrayfun(@(x) x.ureventinfo.latency(1), data);
                                else
                                    aas_log(aap,true,'ERROR: original eventinfo (ureventinfo) is not available, use EEGLAB dataset as input')
                                end
                                dat = cellfun(@(x) x.powspctrm, tf, 'UniformOutput', false);
                                timefreqModel.powspctrm = cat(3,dat{:});
                                dat = cellfun(@(x) x.crsspctrm, tf, 'UniformOutput', false);
                                timefreqModel.crsspctrm = cat(3,dat{:});
                        end
                    else
                        clear tf
                        for i = 1:numel(data)
                            cfg = tfacfg;
                            cfg.trials = find(data(i).trialinfo==trialinfo);
                            cfg.keeptrials = 'yes';
                            tmptfr = ft_freqanalysis(cfg, data(i));
                            X = m.samplevector(data(i).sampleinfo(:,1)); X = X - mean(X); X = X / (max(X)-min(X));
                            X = [X ones(size(X,1),1)];
                            for ch = 1:size(tmptfr.powspctrm,2)
                                for f = 1:size(tmptfr.powspctrm,3)
                                    for t = 1:size(tmptfr.powspctrm,4)
                                        b = X\tmptfr.powspctrm(:,ch,f,t);
                                        resavg(ch,f,t) = b(1);
                                        resvar(ch,f,t) = var(X*b - tmptfr.powspctrm(:,ch,f,t));
                                        resdof(ch,f,t) = size(tmptfr.powspctrm,1) - size(X,2);
                                    end
                                end
                            end
                            tf{i} = timefreqModel;
                            tf{i}.powspctrm = resavg;
                            tf{i}.var = resvar;
                            tf{i}.dof = resdof;
                        end
                        cfg = combinecfg;
                        cfg.normalise = 'yes';
                        cfg.weights = weights;
                        timefreqModel = ft_combine(cfg,tf{:});
                    end
                          
                    meeg_diagnostics_TFR(timefreqModel,diag,[m.name '_' eventLabel],fullfile(aas_getsesspath(aap,subj,sess),['diagnostic_' mfilename  '_' m.name '_' eventLabel]));

                    conEvents{e} = timefreqModel;
                    timefreq.(eventLabel).main = timefreqMain;
                    timefreq.(eventLabel).model = timefreqModel;
                end
                if any(cellfun(@(x) isempty(x), conEvents)), continue; end
                
                % contarst events
                cfg = combinecfg; 
                cfg.weights = m.event.weights;
                if prod(cfg.weights) < 0, cfg.contrast = aas_getsetting(aap,'contrastoperation'); end % differential contrast
                timefreq = ft_combine(cfg,conEvents{:});
                meeg_diagnostics_TFR(timefreq,diag,m.name,fullfile(aas_getsesspath(aap,subj,sess),['diagnostic_' mfilename  '_' m.name '_eventcontrast']));
                
                % save/update output
                timefreqFn = fullfile(aas_getsesspath(aap,subj,sess),['timefreq_' m.name '.mat']);
                save(timefreqFn,'timefreq');
                
                % append to stream
                outputFn = {};
                outstreamFn = aas_getoutputstreamfilename(aap,'meeg_session',[subj, sessnum],'timefreq');
                if exist(outstreamFn,'file')
                    outputFn = cellstr(aas_getfiles_bystream(aap,'meeg_session',[subj, sessnum],'timefreq','output'));
                end
                outputFn{end+1} = timefreqFn;
                aap = aas_desc_outputs(aap,'meeg_session',[subj,sess],'timefreq',outputFn);
                
                conSessions{sess} = timefreq;
            end            
            if any(cellfun(@(x) isempty(x), conSessions)), continue; end
            
            % contrast sessions
            cfg = combinecfg; 
            cfg.weights = m.session.weights;
            timefreq = ft_combine(cfg,conSessions{:});
            meeg_diagnostics_TFR(timefreq,diag,m.name,fullfile(aas_getsubjpath(aap,subj),['diagnostic_' mfilename  '_' m.name]));
            
            % save/update output
            timefreq.cfg = []; % remove provenance to save space
            timefreqFn = fullfile(aas_getsubjpath(aap,subj),['timefreq_' m.name '.mat']);
            save(timefreqFn,'timefreq');
            
            % append to stream
            outputFn = {};
            outstreamFn = aas_getoutputstreamfilename(aap,'subject',subj,'timefreq');
            if exist(outstreamFn,'file')
                outputFn = cellstr(aas_getfiles_bystream(aap,'subject',subj,'timefreq','output'));
            end
            outputFn{end+1} = timefreqFn;
            aap = aas_desc_outputs(aap,'subject',subj,'timefreq',outputFn);
        end
        
        FT.rmExternal('spm12');
        FT.unload;
    case 'checkrequirements'
        if ~aas_cache_get(aap,'eeglab'), aas_log(aap,false,'EEGLAB is not found -> You will not be able process EEGLAB data'); end
        if ~aas_cache_get(aap,'fieldtrip'), aas_log(aap,true,'FieldTrip is not found'); end
end
end