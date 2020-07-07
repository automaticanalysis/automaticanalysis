function [aap, resp] = aamod_meeg_timelockanalysis(aap,task,subj)

resp='';

switch task
    case 'report'
        models = strrep(spm_file(cellstr(aas_getfiles_bystream(aap,'subject',subj,'timelock')),'basename'),'timelock_','')';
        aap = aas_report_add(aap,subj,'<table><tr>');
        for m = models, aap = aas_report_add(aap,subj,['<th>Model: ' m{1} '</th>']); end
        aap = aas_report_add(aap,subj,'</tr><tr>');
        
        for m = models
            aap = aas_report_add(aap,subj,'<td>');
            
            aap=aas_report_addimage(aap,subj,fullfile(aas_getsubjpath(aap,subj),['diagnostic_' mfilename  '_' m{1} '_multiplot.jpg']));
            aap=aas_report_addimage(aap,subj,fullfile(aas_getsubjpath(aap,subj),['diagnostic_' mfilename  '_' m{1} '_topoplot.jpg']));
            aap=aas_report_addimage(aap,subj,fullfile(aas_getsubjpath(aap,subj),['diagnostic_' mfilename  '_' m{1} '_topoplot.avi']));
            
            aap = aas_report_add(aap,subj,'</td>');
        end
        
        aap = aas_report_add(aap,subj,'</tr></table>');
    case 'doit'
        [junk, FT] = aas_cache_get(aap,'fieldtrip');
        FT.load;
        FT.addExternal('spm12');
        
        diag = aas_getsetting(aap,'diagnostics');
        
        models = aas_getsetting(aap,'trialmodel');
        subjmatches=strcmp(aap.acq_details.subjects(subj).subjname,{models.subject});
        if ~any(subjmatches), aas_log(aap,true,['No trialmodel specification found for ' aas_getsubjdesc(aap,subj)]); end
        peakdefs = aas_getsetting(aap,'peakdef');
        peakdefs = peakdefs(strcmp(aap.acq_details.subjects(subj).subjname,{peakdefs.subject}));
        
        combinecfg = [];
        combinecfg.parameter = 'avg';
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
                            data(seg) = dat.data;
                        case 'eeglab'
                            FT.unload;
                            if seg == 1
                                [junk, EL] = aas_cache_get(aap,'eeglab');
                                EL.load;
                            else
                                EL.reload;
                            end
                            EEG = pop_loadset(meegfn{seg});
                            EL.unload;
                            FT.reload;                            
                            data(seg) = eeglab2fieldtripER(EEG);
                    end
                end
                if isfield(data,'ursamplenum'), data = rmfield(data,'ursamplenum'); end
                
                % process events
                kvs = regexp(spm_file(meegfn{1},'basename'),'[A-Z]+-[0-9]+','match');
                events = cellfun(@(x) strsplit(x,'-'), kvs,'UniformOutput',false);
                conEvents = cell(1,numel(m.event.names));
                
                for e = 1:numel(m.event.names)
                    eventLabel = m.event.names{e};
                    trialinfo = str2double(events{cellfun(@(x) strcmp(x{1},eventLabel), events)}{2});
                    aas_log(aap,false,sprintf('INFO: processing event %s with trialinfo %d',eventLabel,trialinfo));
                    
                    % main 
                    tl = {}; weights = [];
                    for i = 1:numel(data)
                        cfg = [];
                        cfg.trials = find(data(i).trialinfo==trialinfo);
                        if isempty(cfg.trials), continue; end
                        tl{end+1} = ft_timelockanalysis(cfg, data(i));
                        if aas_getsetting(aap,'weightedaveraging')
                            weights(end+1) = numel(cfg.trials);
                        else
                            weights(end+1) = 1;
                        end
                    end
                    if isempty(tl), continue; end
                    
                    cfg = combinecfg; 
                    cfg.normalise = 'yes';
                    cfg.weights = weights;
                    timelockMain = ft_combine(cfg,tl{:});
                    
                    diagFn = fullfile(aas_getsesspath(aap,subj,sess),['diagnostic_' mfilename '_' eventLabel]);
                    if ~exist([diagFn '_multiplot.jpg'],'file')
                        meeg_diagnostics_ER(timelockMain,diag,eventLabel,diagFn);
                    end
                    
                    % trialmodel
                    timelockModel = timelockMain;
                    if ischar(m.samplevector) && strcmp(m.samplevector,'avg') % average - it is done
                        % do nothing
                    else
                        clear tl
                        for i = 1:numel(data)
                            cfg = [];
                            cfg.trials = find(data(i).trialinfo==trialinfo);
                            cfg.keeptrials = 'yes';
                            tmper = ft_timelockanalysis(cfg, data(i));
                            X = m.samplevector(tmper.sampleinfo(:,1)); X = X - mean(X); X = X / (max(X)-min(X));
                            X = [X ones(size(X,1),1)];
                            for ch = 1:size(tmper.trial,2)
                                for t = 1:size(tmper.trial,3)
                                    b = X\tmper.trial(:,ch,t);
                                    resavg(ch,t) = b(1);
                                    resvar(ch,t) = var(X*b - tmper.trial(:,ch,t));
                                    resdof(ch,t) = size(tmper.trial,1) - size(X,2);
                                end
                            end
                            tl{i} = timelockModel;
                            tl{i}.avg = resavg;
                            tl{i}.var = resvar;
                            tl{i}.dof = resdof;
                        end
                        cfg = combinecfg;
                        cfg.normalise = 'yes';
                        cfg.weights = weights;
                        timelockModel = ft_combine(cfg,tl{:});
                    end
                          
                    meeg_diagnostics_ER(timelockModel,diag,[m.name '_' eventLabel],fullfile(aas_getsesspath(aap,subj,sess),['diagnostic_' mfilename  '_' m.name '_' eventLabel]));

                    conEvents{e} = timelockModel;
                    timelock.(eventLabel).main = timelockMain;
                    timelock.(eventLabel).model = timelockModel;
                end
                if any(cellfun(@(x) isempty(x), conEvents)), continue; end
                
                % contarst events
                cfg = combinecfg; 
                cfg.weights = m.event.weights;
                timelockCon = ft_combine(cfg,conEvents{:});
                meeg_diagnostics_ER(timelockModel,diag,m.name,fullfile(aas_getsesspath(aap,subj,sess),['diagnostic_' mfilename  '_' m.name '_eventcontrast']));
  
                timelock.contrast = timelockCon;
                
                % save/update output
                timelockFn = fullfile(aas_getsesspath(aap,subj,sess),['timelock_' m.name '.mat']);
                save(timelockFn,'timelock');
                
                % append to stream
                outputFn = {};
                outstreamFn = aas_getoutputstreamfilename(aap,'meeg_session',[subj, sessnum],'timelock');
                if exist(outstreamFn,'file')
                    outputFn = cellstr(aas_getfiles_bystream(aap,'meeg_session',[subj, sessnum],'timelock','output'));
                end
                outputFn{end+1} = timelockFn;
                aap = aas_desc_outputs(aap,'meeg_session',[subj,sess],'timelock',outputFn);
                
                conSessions{sess} = timelock.contrast;                 
            end            
            if any(cellfun(@(x) isempty(x), conSessions)), continue; end
            
            % contrast sessions
            cfg = combinecfg; 
            cfg.weights = m.session.weights;
            timelock = ft_combine(cfg,conSessions{:});
            meeg_diagnostics_ER(timelockModel,diag,m.name,fullfile(aas_getsubjpath(aap,subj),['diagnostic_' mfilename  '_' m.name]));
            
            % save/update output
            timelockFn = fullfile(aas_getsubjpath(aap,subj),['timelock_' m.name '.mat']);
            save(timelockFn,'timelock');
            
            % append to stream
            outputFn = {};
            outstreamFn = aas_getoutputstreamfilename(aap,'subject',subj,'timelock');
            if exist(outstreamFn,'file')
                outputFn = cellstr(aas_getfiles_bystream(aap,'subject',subj,'timelock','output'));
            end
            outputFn{end+1} = timelockFn;
            aap = aas_desc_outputs(aap,'subject',subj,'timelock',outputFn);
            
            % detect peak
            if ~isempty(peakdefs)
                trialmatches = arrayfun(@(x) any(strcmp(x.trial,m.name)), peakdefs);
                if any(trialmatches)
                    for p = peakdefs(trialmatches).peakdef
                        cfg = [];
                        cfg.inflection = p.direction;
                        cfg.latency = p.toi/1000;
                        cfg.neighbourwidth = aas_getsetting(aap,'peakneighbours')/1000;
                        peak = ft_detectpeak(cfg,timelock);
                        
                        peakFn = fullfile(aas_getsubjpath(aap,subj),['peak_' m.name '_' p.name '.mat']);
                        save(peakFn,'peak');
                        
                        % append to stream
                        outputFn = {};
                        outstreamFn = aas_getoutputstreamfilename(aap,'subject',subj,'peak');
                        if exist(outstreamFn,'file')
                            outputFn = cellstr(aas_getfiles_bystream(aap,'subject',subj,'peak','output'));
                        end
                        outputFn{end+1} = peakFn;
                        aap = aas_desc_outputs(aap,'subject',subj,'peak',outputFn);
                    end
                end
            end
        end
        
        FT.rmExternal('spm12');
        FT.unload;
    case 'checkrequirements'
        if ~aas_cache_get(aap,'eeglab'), aas_log(aap,false,'EEGLAB is not found -> You will not be able process EEGLAB data'); end
        if ~aas_cache_get(aap,'fieldtrip'), aas_log(aap,true,'FieldTrip is not found'); end
end
end