
function [aap, resp] = aamod_meeg_statistics(aap,task)

resp='';

switch task
    case 'report'
        RES = {'multiplot.jpg', 'topoplot.jpg', 'topoplot.avi'};
        
        aap = aas_report_add(aap,[],'<h3>Electrode neighbourhood</h3>');
        aap = aas_report_addimage(aap,[],fullfile(aas_getstudypath(aap),['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_neighbours.jpg']));
        
        models = aas_getsetting(aap,'model'); models(1) = []; 
        
        aap = aas_report_add(aap,[],'<table><tr>');
        for m = models, aap = aas_report_add(aap,[],['<th>Model: ' m.name '</th>']); end
        aap = aas_report_add(aap,[],'</tr><tr>');
        
        for m = models
            aap = aas_report_add(aap,[],'<td valign="top">');
            
            savepath = fullfile(aas_getstudypath(aap),['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_' m.name]);
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
               
        inpstreams = aas_getstreams(aap,'input');
        for subj = 1:numel(aap.acq_details.subjects)
            allFnTL{subj} = cellstr(aas_getfiles_bystream(aap,'subject',subj,inpstreams{1}));
            if aas_isfile_bystream(aap,'subject',subj,'peak')
                allFnP{subj} = cellstr(aas_getfiles_bystream(aap,'subject',subj,'peak'));
            end
        end
        
        % define neighbours based on first available data
        cfg = [];
        cfg.method      = 'triangulation'; 
        cfg.feedback    = 'yes';
        dat = load(allFnTL{1}{1});
        neighbours = ft_prepare_neighbours(cfg, dat.(char(fieldnames(dat))));
        set(gcf,'position',[0,0,1080 1080]);
        set(gcf,'PaperPositionMode','auto');
        print(gcf,'-noui',fullfile(aas_getstudypath(aap),['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_neighbours']),'-djpeg','-r300');
        close(gcf);
        
        thr = aas_getsetting(aap,'threshold');
        statcfg = [];
        statcfg.channel   = 'all';
        statcfg.avgovertime = 'no';
        switch inpstreams{1}
            case 'timelock'
                statcfg.parameter   = 'avg';
                fstat = @ft_timelockstatistics;
                fgrandavg = @ft_timelockgrandaverage;
                fdiag = @meeg_diagnostics_ER;
                avglat = 'latency';
            case 'timefreq'
                statcfg.parameter   = 'powspctrm';
                fstat = @ft_freqstatistics;
                fgrandavg = @ft_freqgrandaverage;
                fdiag = @meeg_diagnostics_TFR;
                avglat = 'toilim';
        end
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

        avgcfg = [];
        avgcfg.channel   = 'all';
        avgcfg.(avglat)   = 'all';
                
        statplotcfg = aas_getsetting(aap,'diagnostics');
        statplotcfg.layout = ft_prepare_layout([], dat.(char(fieldnames(dat))));
        
        models = aas_getsetting(aap,'model'); models(1) = []; 
        for m = models
            savepath{1} = fullfile(aas_getstudypath(aap),['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_' m.name]);
            if ~ischar(m.timewindow), m.timewindow = m.timewindow / 1000; end % in seconds
            
            cfg = [];
            cfg{1} = statcfg;
            pcfg = statplotcfg;
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
                        
            allInp = {}; subjmodel = [];
            for subj = 1:numel(m.subjects)
                if ~any(m.trialmodel{1}=='_') % ER/TFR
                    inpType = inpstreams{1};
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
                        dat = dat.(inpType);
                        if isfield(dat,'time')
                            switch aas_getsetting(aap,'selectoverlappingdata.time')
                                case 'auto'
                                    roundtime = mode(diff(dat.time));
                                case 'ignore'
                                    dat.time = 0:numel(dat.time)-1;
                                    roundtime = 1;
                                otherwise
                                    roundtime = aas_getsetting(aap,'selectoverlappingdata.time');
                            end
                            dat.time = round(dat.time/roundtime)*roundtime; % round timing
                        end
                        allInp{end+1} = dat;
                        subjmodel(end+1) = subj;
                    end
                end
            end
            
            if isfield(dat,'time')
                % select only overlapping timepoints
                % - identify overlapping latencies
                allLatencies = cellfun(@(x) x.time, allInp,'UniformOutput',false);
                allLatencies = unique(horzcat(allLatencies{:}));
                allTime = zeros(1,numel(allLatencies));
                for i = 1:numel(allInp)
                    allTime(i,arrayfun(@(x) find(allLatencies==x), allInp{i}.time)) = 1;
                end
                
                if ischar(aas_getsetting(aap,'selectoverlappingdata.subjects')) || aas_getsetting(aap,'selectoverlappingdata.subjects') < 1
                    if strcmp(aas_getsetting(aap,'selectoverlappingdata.time'),'ignore')
                        [ntrial,indsubj] = sort(sum(allTime,2),'descend'); % number of trials as a function of number of subjects
                        if ischar(aas_getsetting(aap,'selectoverlappingdata.subjects')) % auto
                            atrial = ntrial'.*(1:numel(ntrial));
                            i = find(atrial==max(atrial));
                        else
                            for i = ceil(numel(allInp)*aas_getsetting(aap,'selectoverlappingdata.subjects')):numel(allInp)-1 % check if we can add more subjects without decreasing the number of trials
                                if ntrial(i+1) < ntrial(i), break; end
                            end
                        end
                    else
                        aas_log(aap,false,'WARNING: subject selection works only when latencies are ignored')
                        i = numel(allInp);
                    end
                    aas_log(aap,false,sprintf('INFO: %d subject(s) selected\n\tThe final subject list is saved in groupStat{1}.stat.subjects.',i));
                    allTime = allTime(indsubj(1:i),:);
                    allInp = allInp(indsubj(1:i));
                    m.groupmodel = m.groupmodel(indsubj(1:i));
                    subjmodel = subjmodel(indsubj(1:i));
                end
                allTime = find(sum(allTime)==size(allTime,1));
                
                % - select timepoints
                dimTime = find(strcmp(strsplit(allInp{1}.dimord,'_'),'time'));
                if ndims(allInp{i}.(statcfg.parameter)) ~= dimTime
                    ass_log(aap,true,'time is assumed to be the last dimension')
                end

                for i = 1:numel(allInp)
                    allInp{i}.time = allInp{i}.time(allTime);
                    switch dimTime
                        case 1
                            allInp{i}.(statcfg.parameter) = allInp{i}.(statcfg.parameter)(allTime);
                        case 2
                            allInp{i}.(statcfg.parameter) = allInp{i}.(statcfg.parameter)(:,allTime);
                        case 3
                            allInp{i}.(statcfg.parameter) = allInp{i}.(statcfg.parameter)(:,:,allTime);
                        otherwise
                            ass_log(aap,true,'More than three dimensions are not supported')
                    end
                    
                    if isfield(allInp{i},'trialinfo'), allInp{i}.trialinfo = allInp{i}.trialinfo(allTime); end
                end
                
                if strcmp(aas_getsetting(aap,'selectoverlappingdata.time'),'ignore') && ~isempty(pcfg.snapshottwoi) % equally divide trials
                    step = (max(allInp{i}.time) - min(allInp{i}.time))/pcfg.snapshottwoi;
                    pcfg.snapshottwoi = [min(allInp{i}.time):step:(max(allInp{i}.time)-step)]'*1000;
                    pcfg.snapshottwoi(:,2) = pcfg.snapshottwoi(:,1)+step*1000;
                end
            end
            
            cfg{1}.design(1,:) = m.groupmodel;
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
                                    newallInp{i} = ft_math(struct('parameter',cfg{1}.parameter,'operation','subtract'),allInp{origInd});
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
                otherwise % regression
                    aas_log(aap,false,'INFO: regression design detected')
                    if numel(unique(subjmodel)) == 1
                        cfg{1}.statistic   = 'ft_statfun_depsamplesregrT';
                        cfg{1}.uvar = 2;
                        cfg{1}.design(2,:) = 1;
                    elseif numel(unique(subjmodel)) == numel(m.groupmodel)
                        cfg{1}.statistic   = 'ft_statfun_indepsamplesregrT';
                    else
                        aas_log(aap,false,'mixed regression models are not yet implemented')
                        continue;
                    end
            end
            if ~isfield(dat,'time')
                cfg{1}.latency = 'all';
                cfg{1}.correctm = thr.correctiontimepoint;
            end
            if strcmp(inpType, 'peak')  % custom analysis and plotting
                cfg{1}.parameter = 'amp';
                savepath{2} = savepath{1};
                savepath{1} = [savepath{1} '_amp'];
                cfg{2} = cfg{1};
                cfg{2}.parameter = 'lat';
                savepath{2} = [savepath{2} '_lat'];
            end
            
            for c = 1:numel(cfg)
                statFn = [savepath{c} '_' cfg{c}.correctm];
                
                stat = fstat(cfg{c}, allInp{:});
                
                % plot
                groupStat = {};
                avgcfg.parameter = cfg{c}.parameter;
                avgcfg.(avglat) = cfg{c}.latency;
                if isfield(pcfg,'snapshottwoi') && ~isempty(pcfg.snapshottwoi)
                    pcfg.snapshottwoi = pcfg.snapshottwoi(...
                        pcfg.snapshottwoi(:,1)/1000 >= stat.time(1) & ...
                        pcfg.snapshottwoi(:,2)/1000 <= stat.time(end) ...
                        ,:);
                end
                
                if numel(unique(m.groupmodel)) <= 2
                    for g = unique(m.groupmodel)
                        groupStat{end+1} = fgrandavg(avgcfg, allInp{cfg{c}.design(1,:)==g});
                    end
                else
                    groupStat{1} = fgrandavg(avgcfg, allInp{:});
                end
                groupStat{1}.stat = stat;
                fdiag(groupStat,pcfg,m.name,statFn);
                            
                groupStat{1}.stat.subjects = m.subjects(subjmodel);
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