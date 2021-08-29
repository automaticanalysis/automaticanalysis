function [aap, resp] = aamod_meeg_crossfrequencyanalysis(aap,task,subj)
resp='';

switch task
    case 'report'
%         channels = aas_getsetting(aap,'crossfrequencyanalysis.chanphase'); if isempty(channels), channels = {}; end
%         models = strrep(spm_file(cellstr(aas_getfiles_bystream(aap,'subject',subj,'crossfreq')),'basename'),'crossfreq_','')';
%         aap = aas_report_add(aap,subj,'<table id="data"><tr>');
%         aap = aas_report_add(aap,subj,'<th>Phase channel</th>');
%         for m = models, aap = aas_report_add(aap,subj,['<th>Model: ' m{1} '</th>']); end
%         aap = aas_report_add(aap,subj,'</tr>');
% 
%         for ch = channels
%             aap = aas_report_add(aap,subj,'<tr>');
%             aap = aas_report_add(aap,subj,['<td>' ch{1} '</td>']);
%             for m = models
%                 aap = aas_report_add(aap,subj,'<td>');
%                 if isempty(ch{1})
%                     fnimg = fullfile(aas_getsubjpath(aap,subj),['diagnostic_' mfilename  '_' m{1} '_multiplot.jpg']);
%                 else
%                     fnimg = fullfile(aas_getsubjpath(aap,subj),['diagnostic_' mfilename  '_' m{1} '_' ch{1} '_multiplot.jpg']);
%                 end
%                 if exist(fnimg,'file'), aap=aas_report_addimage(aap,subj,fnimg); end
%                 aap = aas_report_add(aap,subj,'</td>');
%             end
%             aap = aas_report_add(aap,subj,'</tr>');
%         end
% 
%         aap = aas_report_add(aap,subj,'</table>');
    case 'doit'
        [~, FT] = aas_cache_get(aap,'fieldtrip');
        FT.load;
        FT.addExternal('spm12');

        cfa = aas_getsetting(aap,'crossfrequencyanalysis');
        cfacfg = keepfields(cfa,'method');
        if ~isempty(cfa.chanphase), cfacfg.chanlow = cfa.chanphase; end
        if ~isempty(cfa.chanamp), cfacfg.chanhigh = cfa.chanamp; end
        cfacfg.freqlow = cfa.foiphase([1 end]);
        cfacfg.freqhigh = cfa.foiamp([1 end]);
        if ~isempty(cfa.nphasebins), cfacfg.nphase = cfa.nphasebins; end
        cfacfg.keeptrials  = 'no';
               
        tfacfg.method	 = 'mtmconvol';
        tfacfg.precision = 'single';
        tfacfg.pad       = 'nextpow2';
        tfacfg.output    = 'fourier';
        tfacfg.taper     = 'hanning';
        tfacfg.toi       = 'all';
        tfacfg.keeptrials  = 'yes';
        tfacfg.foi = [cfa.foiphase cfa.foiamp];
        tfacfg.t_ftimwin = aas_getsetting(aap,'timefrequencyanalysis.twoicps')./tfacfg.foi;
        
        combinecfg = [];
        combinecfg.parameter = 'crsspctrm';
        combinecfg.normalise = 'no';
        
        diagcfg = aas_getsetting(aap,'diagnostics');
        diagcfg.parameter = 'crsspctrm';
              
        models = aas_getsetting(aap,'trialmodel');
        subjmatches=strcmp(aap.acq_details.subjects(subj).subjname,{models.subject});
        if ~any(subjmatches), aas_log(aap,true,['No trialmodel specification found for ' aas_getsubjdesc(aap,subj)]); end
        
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
                            f = fieldnames(dat);
                            data(seg) = ft_struct2single(dat.(f{1}));
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
                            FT.reload;                            
                            data(seg) = ft_struct2single(eeglab2fieldtripER(EEG,'reorient',1));
                            EL.unload; % FieldTrip's eeglab toolbox is incomplete
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
                    cf = {}; weights = [];
                    for i = 1:numel(data)
                        cfg = tfacfg;
                        cfg.trials = find(data(i).trialinfo==trialinfo);
                        if isempty(cfg.trials), continue; end
                        tf = ft_freqanalysis(cfg, data(i));
                        tf.dimord = strrep(tf.dimord,'rpttap','rpt'); % ft_selectdata does not recognise rpttap
                        cfacfg.freqlow = arrayfun(@(x) tf.freq(abs(tf.freq-x)==min(abs(tf.freq-x))),cfacfg.freqlow);
                        cfacfg.freqhigh = arrayfun(@(x) tf.freq(abs(tf.freq-x)==min(abs(tf.freq-x))),cfacfg.freqhigh);
                        tmpcf = ft_crossfrequencyanalysis(cfacfg,tf);
                        indnoCF = arrayfun(@(x) all(~nanfill(squeeze(tmpcf.crsspctrm(x,:,:,:)),0),'all'),1:size(tmpcf.labelcmb,1));
                        tmpcf.crsspctrm(indnoCF,:,:,:) = [];
                        if isfield(tmpcf,'labelcmb'), tmpcf.labelcmb(indnoCF,:) = [];
                        else, tmpcf.label(indnoCF,:) = []; end
                        cf{end+1} = tmpcf;
                        cf{end}.elec = tf.elec;
                        if aas_getsetting(aap,'weightedaveraging')
                            weights(end+1) = numel(cfg.trials);
                        else
                            weights(end+1) = 1;
                        end
                    end
                    if isempty(cf), continue; end
                    
                    % combine
                    cfg = combinecfg; 
                    cfg.normalise = 'yes';
                    cfg.weights = weights;
                    crossfreqMain = ft_combine(cfg,cf{:});
                    
                    meeg_diagnostics_CF(crossfreqMain,diagcfg,eventLabel,fullfile(aas_getsesspath(aap,subj,sess),['diagnostic_' mfilename  '_' eventLabel]));
                    
                    % trialmodel
                    crossfreqModel = crossfreqMain;
                    if ischar(m.samplevector)
                        switch m.samplevector
                            case 'avg' % average - it is done
                            % do nothing
                            case 'segmentavg'
                                crossfreqModel.dimord    = [crossfreqModel.dimord '_time'];
                                if isfield(data(1),'ureventinfo')
                                    crossfreqModel.time = arrayfun(@(x) x.ureventinfo.latency(1), data);
                                else
                                    aas_log(aap,true,'ERROR: original eventinfo (ureventinfo) is not available, use EEGLAB dataset as input')
                                end
                                dat = cellfun(@(x) x.crsspctrm, cf, 'UniformOutput', false);
                                crossfreqModel.crsspctrm = cat(ndims(dat{1})+1,dat{:});
                        end
                    else
                        clear cf
                        for i = 1:numel(data)
                            cfg = tfacfg;
                            cfg.trials = find(data(i).trialinfo==trialinfo);
                            tf = ft_freqanalysis(cfg, data(i));
                            tf.dimord = strrep(tf.dimord,'rpttap','rpt'); % ft_selectdata does not recognise rpttap
                            cfg = cfacfg;
                            cfg.keeptrials = 'yes';
                            cfg.freqlow = arrayfun(@(x) tf.freq(abs(tf.freq-x)==min(abs(tf.freq-x))),cfg.freqlow);
                            cfg.freqhigh = arrayfun(@(x) tf.freq(abs(tf.freq-x)==min(abs(tf.freq-x))),cfg.freqhigh);
                            tmpcf = ft_crossfrequencyanalysis(cfg,tf);                                                    
                            indnoCF = arrayfun(@(x) all(~nanfill(squeeze(tmpcf.crsspctrm(x,:,:,:)),0),'all'),1:size(tmpcf.labelcmb,1));
                            tmpcf.crsspctrm(indnoCF,:,:,:,:) = [];
                            if isfield(tmpcf,'labelcmb'), tmpcf.labelcmb(indnoCF,:) = [];
                            else, tmpcf.label(indnoCF,:) = []; end
                            tmpcf.elec = tf.elec;
                            
                            X = m.samplevector(data(i).sampleinfo(:,1)); X = X - mean(X); X = X / (max(X)-min(X));
                            X = [X ones(size(X,1),1)];
                            for ch = 1:size(tmptfr.crsspctrm,2)
                                for f1 = 1:size(tmptfr.crsspctrm,3)
                                    for f2 = 1:size(tmptfr.crsspctrm,4)
                                        for ph = 1:size(tmptfr.crsspctrm,5)
                                            b = X\tmptfr.crsspctrm(:,ch,f1,f2,ph);
                                            resavg(ch,f1,f2,ph) = b(1);
                                            resvar(ch,f1,f2,ph) = var(X*b - tmptfr.powspctrm(:,ch,f1,f2));
                                            resdof(ch,f1,f2,ph) = size(tmptfr.powspctrm,1) - size(X,2);
                                        end
                                    end
                                end
                            end
                            cf{i} = crossfreqModel;
                            cf{i}.powspctrm = resavg;
                            cf{i}.var = resvar;
                            cf{i}.dof = resdof;
                        end
                        
                        % combine
                        cfg = combinecfg;
                        cfg.normalise = 'yes';
                        cfg.weights = weights;
                        crossfreqModel = ft_combine(cfg,cf{:});
                    end
                          
                    meeg_diagnostics_CF(crossfreqModel,diagcfg,[m.name '_' eventLabel],fullfile(aas_getsesspath(aap,subj,sess),['diagnostic_' mfilename  '_' m.name '_' eventLabel]));

                    conEvents{e} = crossfreqModel;
                    crossfreq.(eventLabel).main = crossfreqMain;
                    crossfreq.(eventLabel).model = crossfreqModel;
                end
                if any(cellfun(@(x) isempty(x), conEvents)), continue; end
                
                % contarst events
                cfg = combinecfg; 
                cfg.weights = m.event.weights;
                if prod(weights) < 0, cfg.contrast = aas_getsetting(aap,'contrastoperation'); end % differential contrast
                crossfreq = ft_combine(cfg,conEvents{:});
                meeg_diagnostics_CF(crossfreq,diagcfg,m.name,fullfile(aas_getsesspath(aap,subj,sess),['diagnostic_' mfilename  '_' m.name '_eventcontrast']));
              
                % save/update output
                crossfreqFn = fullfile(aas_getsesspath(aap,subj,sess),['crossfreq_' m.name '.mat']);
                save(crossfreqFn,'crossfreq');
                
                % append to stream
                outputFn = {};
                outstreamFn = aas_getoutputstreamfilename(aap,'meeg_session',[subj, sessnum],'crossfreq');
                if exist(outstreamFn,'file')
                    outputFn = cellstr(aas_getfiles_bystream(aap,'meeg_session',[subj, sessnum],'crossfreq','output'));
                end
                outputFn{end+1} = crossfreqFn;
                aap = aas_desc_outputs(aap,'meeg_session',[subj,sess],'crossfreq',outputFn);
                
                conSessions{sess} = crossfreq;
            end            
            if any(cellfun(@(x) isempty(x), conSessions)), continue; end
            
            % contrast sessions
            cfg = combinecfg; 
            cfg.weights = m.session.weights;
            crossfreq = ft_combine(cfg,conSessions{:});
            meeg_diagnostics_CF(crossfreq,diagcfg,m.name,fullfile(aas_getsubjpath(aap,subj),['diagnostic_' mfilename  '_' m.name]));
           
            % save/update output
            crossfreq.cfg = []; % remove provenance to save space
            crossfreqFn = fullfile(aas_getsubjpath(aap,subj),['crossfreq_' m.name '.mat']);
            save(crossfreqFn,'crossfreq');
            
            % append to stream
            outputFn = {};
            outstreamFn = aas_getoutputstreamfilename(aap,'subject',subj,'crossfreq');
            if exist(outstreamFn,'file')
                outputFn = cellstr(aas_getfiles_bystream(aap,'subject',subj,'crossfreq','output'));
            end
            outputFn{end+1} = crossfreqFn;
            aap = aas_desc_outputs(aap,'subject',subj,'crossfreq',outputFn);
        end
        
        FT.rmExternal('spm12');
        FT.unload;
    case 'checkrequirements'
        if ~aas_cache_get(aap,'eeglab'), aas_log(aap,false,'EEGLAB is not found -> You will not be able process EEGLAB data'); end
        if ~aas_cache_get(aap,'fieldtrip'), aas_log(aap,true,'FieldTrip is not found'); end
        
        cfa = aas_getsetting(aap,'crossfrequencyanalysis');
        if strcmp(cfa.method, 'pac') && isempty(cfa.nphasebins)
            aas_log(aap,true,'ERROR: PAC requires nphasebins');
        end
end
end

function t = nanfill(dat,v)
t = dat;
t(isnan(t)) = v;
end