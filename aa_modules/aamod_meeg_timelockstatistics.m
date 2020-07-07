function [aap, resp] = aamod_meeg_timelockstatistics(aap,task)

resp='';

switch task
    case 'report'
        RES = {'multiplot.jpg', 'topoplot.jpg', 'topoplot.avi'};
        
        aap = aas_report_add(aap,[],'<h3>Electrode neighbourhood</h3>');
        aap = aas_report_addimage(aap,[],fullfile(aas_getstudypath(aap),['diagnostic_' mfilename '_neighbours.jpg']));
        
        models = aas_getsetting(aap,'model'); models(1) = []; 
        
        aap = aas_report_add(aap,[],'<table><tr>');
        for m = models, aap = aas_report_add(aap,[],['<th>Model: ' m.name '</th>']); end
        aap = aas_report_add(aap,[],'</tr><tr>');
        
        for m = models
            aap = aas_report_add(aap,[],'<td valign="top">');
            
            savepath = fullfile(aas_getstudypath(aap),['diagnostic_' mfilename '_' m.name]);
            res = cellstr(spm_select('FPList',spm_file(savepath,'path'),['^' spm_file(savepath,'basename') '_.*']));
            for r = RES
                ind = cell_index(res,r{1});
                if ind
                    for i = ind % peaks have multiple hit (amp+lat)
                        aap = aas_report_addimage(aap,[],res{i});
                    end
                end
            end
                        
            aap = aas_report_add(aap,[],'</td>');
        end
        aap = aas_report_add(aap,[],'</tr></table>');
    case 'doit'
        [junk, FT] = aas_cache_get(aap,'fieldtrip');
        FT.load;
                
        for subj = 1:numel(aap.acq_details.subjects)
            allFnTL{subj} = cellstr(aas_getfiles_bystream(aap,'subject',subj,'timelock'));
            if aas_isfile_bystream(aap,'subject',subj,'peak')
                allFnP{subj} = cellstr(aas_getfiles_bystream(aap,'subject',subj,'peak'));
            end
        end
        
        % define neighbours based on first available data
        cfg = [];
        cfg.method      = 'triangulation'; 
        cfg.feedback    = 'yes';
        dat = load(allFnTL{1}{1});
        neighbours = ft_prepare_neighbours(cfg, dat.timelock);
        set(gcf,'position',[0,0,1080 1080]);
        set(gcf,'PaperPositionMode','auto');
        print(gcf,'-noui',fullfile(aas_getstudypath(aap),['diagnostic_' mfilename '_neighbours']),'-djpeg','-r300');
        close(gcf);
        
        avgcfg = [];
        avgcfg.channel   = 'all';
        avgcfg.latency   = 'all';
        
        thr = aas_getsetting(aap,'threshold');
        statcfg = [];
        statcfg.channel   = 'all';
        statcfg.avgovertime = 'no';
        statcfg.parameter   = 'avg';
        statcfg.alpha       = thr.p;
        statcfg.tail        = 0; % two-tailed
        statcfg.correcttail = 'prob';
        statcfg.method              = thr.method;
        statcfg.correctm            = thr.correctiontimeseries;
        statcfg.clusteralpha        = thr.p;
        statcfg.numrandomization    = thr.iteration;
        statcfg.minnbchan           = thr.neighbours;      % minimal number of neighbouring channels
        statcfg.neighbours          = neighbours; % defined as above
        statcfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable        

        statplotcfg = aas_getsetting(aap,'diagnostics');
        statplotcfg.layout = ft_prepare_layout([], dat.timelock);
        
        models = aas_getsetting(aap,'model'); models(1) = []; 
        for m = models
            savepath{1} = fullfile(aas_getstudypath(aap),['diagnostic_' mfilename '_' m.name]);
            if ~ischar(m.timewindow), m.timewindow = m.timewindow / 1000; end % in seconds
            
            cfg = [];
            cfg{1} = statcfg;
            if ischar(m.channels) || iscell(m.channels)
                cfg{1}.channel = m.channels; % 'all', selected channels, or a single channel 
            elseif isstruct(m.channels) % modelling channel effect
                cfg{1}.channel = {'aa'}; % label for synthetic channel
                combinecfg = [];
                combinecfg.parameter = cfg{1}.parameter;
                combinecfg.weights = m.channels.weights;
                combinecfg.normalise = 'yes';                
            else, aas_log(aap,true,['Wrong channel specification for model ' m.name]);
            end
            cfg{1}.latency = m.timewindow;
            cfg{1}.design(1,:) = m.groupmodel;
            
            allInp = {}; subjmodel = [];
            for subj = 1:numel(m.subjects)
                if ~any(m.trialmodel{1}=='_') % ER
                    inpType = 'timelock';
                    sFn = allFnTL{cell_index({aap.acq_details.subjects.subjname},m.subjects{subj})};
                else % peak
                    inpType = 'peak';
                    sFn = allFnP{cell_index({aap.acq_details.subjects.subjname},m.subjects{subj})};
                end                
                for iFn = find(sum(cell2mat(cellfun(@(y) cellfun(@(x) ~isempty(x), regexp(spm_file(sFn,'basename'),[inpType '_' y '$'])), m.trialmodel,'UniformOutput',false)),2))
                    if isempty(iFn), aas_log(aap,true,'Trialmodels not found'); end
                    for fn = sFn(iFn)'
                        aas_log(aap,false,sprintf('INFO: load %s',fn{1}));
                        dat = load(fn{1});
                        if strcmp(cfg{1}.channel, 'aa') % split data into channels to be able to model them
                            chdata = {};
                            for ch = m.channels.channels
                                chdata{end+1} = ft_selectdata(struct('channel',ch),dat.(inpType));
                                chdata{end}.label = cfg{1}.channel;
                                chdata{end}.elec.label(strcmp(chdata{end}.elec.label,ch)) = cfg{1}.channel;
                            end
                            dat.(inpType) = ft_combine(combinecfg,chdata{:});
                        end
                        allInp{end+1} = dat.(inpType);
                        subjmodel(end+1) = subj;
                    end
                end
            end
            
            switch numel(unique(m.groupmodel))
                case 1 % group average 
                    aas_log(aap,false,'INFO: onesampleT design detected')
                    cfg{1}.statistic   = 'ft_statfun_onesampleT';
                case 2 % 2 groups
                    if isequal(subjmodel(m.groupmodel == 1),subjmodel(m.groupmodel == 2)) % 2x1 group
                        aas_log(aap,false,'INFO: depsamplesT design detected')
                        cfg{1}.statistic = 'ft_statfun_depsamplesT';
                        cfg{1}.uvar = 2;
                        cfg{1}.design(2,:)  = subjmodel;
                    else % Nx2 groups
                        switch numel(subjmodel(m.groupmodel == 1)) / numel(unique(subjmodel(m.groupmodel == 1)))
                            case 1 % 1x2 groups design
                                aas_log(aap,false,'INFO: indepsamplesT design detected')
                                cfg{1}.statistic = 'ft_statfun_indepsamplesT';
                            case 2 % 2x2 groups design (assume full design)
                                aas_log(aap,false,'INFO: 2x2 mixed design design detected -> calculating difference in within-subject differences')
                                newallInp = {}; newdesign = [];
                                for i = 1:numel(subjmodel)/2
                                    origInd = [(i-1)*2+1 i*2];
                                    if diff(subjmodel(origInd))+diff(m.groupmodel(origInd)) > 0
                                        aas_log(aap,true,'Not a full 2x2 design -> NYI');
                                    end
                                    newallInp{i} = ft_math(struct('parameter','avg','operation','subtract'),allInp{origInd});
                                    newdesign(i) = cfg{1}.design(origInd(1)); % both are supposed to be the same
                                end
                                allInp = newallInp;
                                cfg{1}.design = newdesign;
                                cfg{1}.statistic = 'ft_statfun_indepsamplesT';
                            otherwise
                                aas_log(aap,false,'models with more then 2 (repeated) levels are not yet implemented')
                                continue
                        end
                    end
                otherwise
                    aas_log(aap,false,'models with more then 2 (independent) levels are not yet implemented')
                    continue
            end
            if strcmp(inpType, 'peak') % custom analysis and plotting
                cfg{1}.latency = 'all';
                cfg{1}.correctm = thr.correctiontimepoint;
                cfg{1}.parameter = 'amp';
                savepath{2} = savepath{1};
                savepath{1} = [savepath{1} '_amp'];
                cfg{2} = cfg{1};
                cfg{2}.parameter = 'lat';
                savepath{2} = [savepath{2} '_lat'];
            end
            
            for c = 1:numel(cfg)
                statFn = [savepath{c} '_' cfg{c}.correctm];
                
                stat = ft_timelockstatistics(cfg{c}, allInp{:});
                
                % plot
                groupStat = {};
                avgcfg.parameter = cfg{c}.parameter;
                avgcfg.latency = cfg{c}.latency;
                if isfield(statplotcfg,'snapshottwoi')
                    statplotcfg.snapshottwoi = statplotcfg.snapshottwoi(...
                        statplotcfg.snapshottwoi(:,1)/1000 >= stat.time(1) & ...
                        statplotcfg.snapshottwoi(:,2)/1000 <= stat.time(end) ...
                        ,:);
                end
                for g = unique(m.groupmodel)
                    groupStat{end+1} = ft_timelockgrandaverage(avgcfg, allInp{cfg{c}.design(1,:)==g});
                end
                groupStat{1}.stat = stat;
                meeg_diagnostics_ER(groupStat,statplotcfg,m.name,statFn);
                            
                statFn = [statFn '.mat'];
                save(statFn,'groupStat');
                
                % append to stream
                outputFn = {};
                outstreamFn = aas_getoutputstreamfilename(aap,'subject',subj,'groupstat');
                if exist(outstreamFn,'file')
                    outputFn = cellstr(aas_getfiles_bystream(aap,'subject',subj,'groupstat','output'));
                end
                outputFn{end+1} = statFn;
                aap = aas_desc_outputs(aap,'study',[],'groupstat',outputFn);
            end
            
        end
        
        FT.unload;

    case 'checkrequirements'
        if ~aas_cache_get(aap,'fieldtrip'), aas_log(aap,true,'FieldTrip is not found'); end
end
end