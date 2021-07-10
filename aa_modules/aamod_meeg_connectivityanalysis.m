function [aap, resp] = aamod_meeg_connectivityanalysis(aap,task,subj)
resp='';

METHOD_MVAR = struct(...
    'LWR',0,...
    'corr_biased',1,...
    'corr_unbiased',6,...
    'VMparcorr_unbiased',2,...
    'VMparcorr_biased',5,...
    'NSparcorr_unbiased',3,...
    'NSparcorr_biased',7,...
    'ARFIT',10,...
    'BURGV',11 ...
);

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

        mvarcfg = [];
        tfacfg = [];
        conncfg = [];
        combinecfg = [];
        
        conncfg = aas_getsetting(aap,'connectivityanalysis');
        % conncfg.keeptrials  = 'no';
        
        if ~isempty(conncfg.foi), tfacfg.foi = conncfg.foi; end
        conncfg = rmfield(conncfg, 'foi');        
        switch conncfg.method
%                 mvarcfg = aas_getsetting(aap,'multivariateautoregressiveanalysis');
%                 if METHOD_MVAR.(mvarcfg.method) == 0
%                     mvarcfg.method = 'bsmart';
%                 else
%                     mvarcfg.mvarmethod = METHOD_MVAR.(mvarcfg.method);
%                     mvarcfg.method = 'biosig';
%                 end
%                 tfacfg.method = 'mvar';
            case {'granger' 'wpli_debiased'}
                tfacfg.method = 'mtmfft';
                tfacfg.taper  = 'dpss';
                tfacfg.tapsmofrq = aas_getsetting(aap,'timefrequencyanalysis.spectralsmoothing');
                tfacfg.pad = 'maxperlen';
                tfacfg.output = 'fourier';
                combinecfg.parameter = [conncfg.method 'spctrm'];
            otherwise
                aas_log(aap,true,'NYI');
        end               
        
        combinecfg.normalise = 'no';
        
        diagcfg = aas_getsetting(aap,'diagnostics');
        diagcfg.parameter = combinecfg.parameter;
              
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
                data = rmfield(data,'included');
                
                % process channels
                if ~isfield(conncfg,'channelcmb')
                    if isempty(conncfg.channels)
                        conncfg.channelcmb = ft_channelcombination({data(1).label,data(1).label},data(1).label,1,2);
                    else
                        if size(conncfg.channels,2) == 1 % channel list
                            channels = intersect(conncfg.channels,data(1).label);
                            conncfg.channelcmb = ft_channelcombination({channels,channels},channels,1,2);
                        else % channel combination
                            channelcmb = conncfg.channels;
                            channelcmb = channelcmb(arrayfun(@(cmb) numel(intersect(channelcmb(cmb,:),data(1).label))==2, 1:size(channelcmb,1)),:);
                            conncfg.channelcmb = channelcmb;
                        end
                    end
                    conncfg.channel = unique(conncfg.channelcmb(:));
                    if ~isempty(mvarcfg)
                        cfg = []; cfg.channel = conncfg.channel;
                        data = arrayfun(@(d) ft_selectdata(cfg,d), data);
                    end
                end
                
                % process events
                kvs = regexp(spm_file(meegfn{1},'basename'),'[A-Z]+-[0-9]+','match');
                events = cellfun(@(x) strsplit(x,'-'), kvs,'UniformOutput',false);
                conEvents = cell(1,numel(m.event.names));
                
                for e = 1:numel(m.event.names)
                    eventLabel = m.event.names{e};
                    trialinfo = str2double(events{cellfun(@(x) strcmp(x{1},eventLabel), events)}{2});
                    aas_log(aap,false,sprintf('INFO: processing event %s with trialinfo %d',eventLabel,trialinfo));
                    
                    conncfg.trials = trialinfo;
                    
                    % main 
                    conn = {}; weights = [];
                    for i = 1:numel(data)                        
                        conn{end+1} = run_connectivity({mvarcfg,tfacfg,conncfg,combinecfg},data(i)); 

                        if aas_getsetting(aap,'weightedaveraging')
                            weights(end+1) = numel(trials);
                        else
                            weights(end+1) = 1;
                        end
                    end
                    indNoConn = cellfun(@(c) isempty(c.(combinecfg.parameter)), conn);
                    conn(indNoConn) = [];
                    weights(indNoConn) = [];
                    
                    if isempty(conn), continue; end
                    
                    % combine
                    cfg = combinecfg; 
                    cfg.normalise = 'yes';
                    cfg.weights = weights;
                    connfreqMain = ft_combine(cfg,conn{:});
                    
                    meeg_diagnostics_conn(connfreqMain,diagcfg,eventLabel,fullfile(aas_getsesspath(aap,subj,sess),['diagnostic_' mfilename  '_' eventLabel]));
                    
                    % trialmodel
                    switch m.samplevector
                        case 'avg' % average - all trials across segments
                            % concatenate data
                            dat = rmfield(data,{'ureventinfo'});
                            dat = num2cell(dat);
                            dat = ft_appenddata(struct('keepsampleinfo','no'),dat{:});
                            
                            connfreqModel = run_connectivity({mvarcfg,tfacfg,conncfg,combinecfg},dat); 
                        case 'segmentavg'
                            connfreqModel = connfreqMain;
                            connfreqModel.dimord    = [connfreqModel.dimord '_time'];
                            if isfield(data(1),'ureventinfo')
                                connfreqModel.time = arrayfun(@(x) x.ureventinfo.latency(1), data);
                                connfreqModel.time(indNoConn) = [];
                            else
                                aas_log(aap,true,'ERROR: original eventinfo (ureventinfo) is not available, use EEGLAB dataset as input')
                            end
                            dat = cellfun(@(x) x.(combinecfg.parameter), conn, 'UniformOutput', false);
                            connfreqModel.(combinecfg.parameter) = cat(ndims(dat{1})+1,dat{:});
                    end
                    
                    meeg_diagnostics_conn(connfreqModel,diagcfg,[m.name '_' eventLabel],fullfile(aas_getsesspath(aap,subj,sess),['diagnostic_' mfilename  '_' m.name '_' eventLabel]));

                    conEvents{e} = connfreqModel;
                    connfreq.(eventLabel).main = connfreqMain;
                    connfreq.(eventLabel).model = connfreqModel;
                end
                if any(cellfun(@(x) isempty(x), conEvents)), continue; end
                
                % contarst events
                cfg = combinecfg; 
                cfg.weights = m.event.weights;
                if prod(weights) < 0, cfg.contrast = aas_getsetting(aap,'contrastoperation'); end % differential contrast
                connfreq = ft_combine(cfg,conEvents{:});
                meeg_diagnostics_conn(connfreq,diagcfg,m.name,fullfile(aas_getsesspath(aap,subj,sess),['diagnostic_' mfilename  '_' m.name '_eventcontrast']));
              
                % save/update output
                connfreqFn = fullfile(aas_getsesspath(aap,subj,sess),['connfreq_' m.name '.mat']);
                save(connfreqFn,'connfreq');
                
                % append to stream
                outputFn = {};
                outstreamFn = aas_getoutputstreamfilename(aap,'meeg_session',[subj, sessnum],'connfreq');
                if exist(outstreamFn,'file')
                    outputFn = cellstr(aas_getfiles_bystream(aap,'meeg_session',[subj, sessnum],'connfreq','output'));
                end
                outputFn{end+1} = connfreqFn;
                aap = aas_desc_outputs(aap,'meeg_session',[subj,sess],'connfreq',outputFn);
                
                conSessions{sess} = connfreq;
            end            
            if any(cellfun(@(x) isempty(x), conSessions)), continue; end
            
            % contrast sessions
            cfg = combinecfg; 
            cfg.weights = m.session.weights;
            connfreq = ft_combine(cfg,conSessions{:});
            meeg_diagnostics_conn(connfreq,diagcfg,m.name,fullfile(aas_getsubjpath(aap,subj),['diagnostic_' mfilename  '_' m.name]));
           
            % save/update output
            connfreq.cfg = []; % remove provenance to save space
            connfreqFn = fullfile(aas_getsubjpath(aap,subj),['connfreq_' m.name '.mat']);
            save(connfreqFn,'connfreq');
            
            % append to stream
            outputFn = {};
            outstreamFn = aas_getoutputstreamfilename(aap,'subject',subj,'connfreq');
            if exist(outstreamFn,'file')
                outputFn = cellstr(aas_getfiles_bystream(aap,'subject',subj,'connfreq','output'));
            end
            outputFn{end+1} = connfreqFn;
            aap = aas_desc_outputs(aap,'subject',subj,'connfreq',outputFn);
        end
        
        FT.rmExternal('spm12');
        FT.unload;
    case 'checkrequirements'
        if ~aas_cache_get(aap,'eeglab'), aas_log(aap,false,'EEGLAB is not found -> You will not be able process EEGLAB data'); end
        if ~aas_cache_get(aap,'fieldtrip'), aas_log(aap,true,'FieldTrip is not found'); end
end
end

function conn = run_connectivity(cfg,dat)
conn = [];
[mvarcfg,tfacfg,conncfg,combinecfg] = cfg{:};

trials = find(dat.trialinfo==conncfg.trials);
conncfg = rmfield(conncfg,'trials');
if isempty(trials), return; end
dat = ft_selectdata(struct('trials',trials),dat);

if ~isempty(mvarcfg), dat = ft_mvaranalysis(mvarcfg,dat); end
if ~isempty(tfacfg), dat = ft_freqanalysis(tfacfg,dat); end
conn = ft_connectivityanalysis(conncfg,dat);

if ~isfield(conn,'labelcmb')
    conn.labelcmb = ft_channelcombination({conn.label,conn.label},conn.label,1,2);
    conn.(combinecfg.parameter) = reshape(conn.(combinecfg.parameter),size(conn.labelcmb,1),[]);
    conn = rmfield(conn,'label');
    conn.dimord = strrep(conn.dimord,'chan_chan','chancmb');
end

indnoCONN = arrayfun(@(x) all(~nanfill(squeeze(conn.(combinecfg.parameter)(x,:,:,:)),0),'all'),1:size(conn.labelcmb,1));
conn.(combinecfg.parameter)(indnoCONN,:,:,:) = [];
conn.labelcmb(indnoCONN,:) = [];
conn.elec = dat.elec;
end

function t = nanfill(dat,v)
t = dat;
t(isnan(t)) = v;
end