
function [aap, resp] = aamod_meeg_statistics(aap,task)

resp='';

switch task
    case 'report'
        RES = {'multiplot.*jpg', 'topoplot.*jpg', 'topoplot.*avi', 'source.*jpg'};
        
        fn = fullfile(aas_getstudypath(aap),['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_neighbours.jpg']);
        if exist(fn,'file')
            aap = aas_report_add(aap,[],'<h3>Electrode neighbourhood</h3>');
            aap = aas_report_addimage(aap,[],fn);
        end
        
        models = aas_getsetting(aap,'model'); models(1) = []; 
        
        aap = aas_report_add(aap,[],'<table><tr>');
        for m = models, aap = aas_report_add(aap,[],['<th>Model: ' m.name '</th>']); end
        aap = aas_report_add(aap,[],'</tr><tr>');
        
        for m = models
            aap = aas_report_add(aap,[],'<td valign="top">');
            
            savepath = fullfile(aas_getstudypath(aap),['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_' m.name]);
            res = cellstr(spm_select('FPList',spm_file(savepath,'path'),['^' spm_file(savepath,'basename') '_.*']));
            for r = RES
                ind = find(cellfun(@(x) ~isempty(regexp(x,r{1},'once')), res))';
                if ~isempty(ind)
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
        dat = load(allFnTL{1}{1}); data = dat.(char(fieldnames(dat)));
        if ~ft_datatype(data,'source')
            neighbours = ft_prepare_neighbours(cfg, data);
            set(gcf,'position',[0,0,720 720]);
            set(gcf,'PaperPositionMode','auto');
            print(gcf,'-noui',fullfile(aas_getstudypath(aap),['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_neighbours']),'-djpeg','-r300');
            close(gcf);
        else
            neighbours = [];
        end
        
        statplotcfg = aas_getsetting(aap,'diagnostics');
        if ~ft_datatype(data,'source')
            statplotcfg.layout = ft_prepare_layout([], dat.(char(fieldnames(dat))));
        end
        
        thr = aas_getsetting(aap,'threshold');
        statcfg = [];
        statcfg.channel   = 'all';
        statcfg.avgovertime = 'no';
        switch inpstreams{1}
            case 'timelock'
                statcfg.parameter   = 'avg';
                fstat = @ft_timelockstatistics;
                fdiag = @meeg_diagnostics_ER;
            case 'timefreq'
                statcfg.parameter   = 'powspctrm';
                statcfg.avgoverfreq = 'yes';
                fstat = @ft_freqstatistics;
                fdiag = @meeg_diagnostics_TFR;
        end
        if ft_datatype(data,'source')
            switch inpstreams{1}
            case 'timefreq'
                statcfg.parameter   = 'pow';
            end
            fstat = @ft_sourcestatistics;
            fdiag = @meeg_diagnostics_source;
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
                        dat.cfg = keepfields(dat.cfg,'included'); % remove provenance to save space
                        allInp{end+1} = ft_struct2single(dat);
                        subjmodel(end+1) = subj;
                    end
                end
            end
            
            if isfield(allInp{1},'time')
                aas_log(aap,false,'INFO: selecting overlapping timepoints');
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
                if isfield(allInp{1},'dimord')
                    dimord = allInp{1}.dimord;
                    par = {statcfg.parameter};
                elseif isfield(allInp{1}.avg,'dimord')
                    dimord = allInp{1}.avg.dimord;
                    par = {'avg' statcfg.parameter};
                else
                    aas_log(aap,true,'no dimension information found');
                end
                dimTime = find(strcmp(strsplit(dimord,'_'),'time'));
                if ndims(getfield(allInp{i}, par{:})) ~= dimTime
                    aas_log(aap,true,'time is assumed to be the last dimension')
                end

                for i = 1:numel(allInp)
                    allInp{i}.time = allInp{i}.time(allTime);
                    dat = getfield(allInp{i}, par{:});
                    switch dimTime
                        case 1
                            allInp{i} = setfield(allInp{i}, par{:}, dat(allTime));
                        case 2
                            allInp{i} = setfield(allInp{i}, par{:}, dat(:,allTime));
                        case 3
                            allInp{i} = setfield(allInp{i}, par{:}, dat(:,:,allTime));
                        otherwise
                            ass_log(aap,true,'More than three dimensions are not supported')
                    end
                    
                    if isfield(allInp{i},'trialinfo'), allInp{i}.trialinfo = allInp{i}.trialinfo(allTime); end
                end
                
                if strcmp(aas_getsetting(aap,'selectoverlappingdata.time'),'ignore') && ~isempty(pcfg.snapshottwoi) % equally divide trials
                    step = (max(allInp{i}.time) - min(allInp{i}.time))/size(pcfg.snapshottwoi,1);
                    pcfg.snapshottwoi = [min(allInp{i}.time):step:(max(allInp{i}.time)-step)]'*1000;
                    pcfg.snapshottwoi(:,2) = pcfg.snapshottwoi(:,1)+step*1000;
                end
            end
            
            if isfield(allInp{1},'pos')
                aas_log(aap,false,'INFO: selecting overlapping positions');
                posJoint = all(cell2mat(cellfun(@(x) reshape(x.cfg.included,[],1), allInp, 'UniformOutput', false)),2);
                allPos = zeros(sum(posJoint),3,numel(allInp));
                for i = 1:numel(allInp)
                    clear dat;
                    dat = allInp{i};
                    toExcl = false(1,size(dat.pos,1));
                    for p = find(reshape(dat.cfg.included,[],1) & ~posJoint)'
                        toExcl(sum(dat.cfg.included(1:p))) = true;
                    end
                    dat.pos(toExcl,:) = [];
                    dat.inside(toExcl,:) = [];
                    dat.avg.(cfg{1}.parameter)(toExcl,:,:) = [];
                    dat.avg.ori(toExcl) = [];
                    if isfield(dat,'tri')
                        [trind, trelem] = ndgrid(1:size(dat.tri,1),1:size(dat.tri,2)); trind = transpose(trind); trelem = transpose(trelem);
                        indtri = logical(sum(arrayfun(@(x,y) any(dat.tri(x,y)==find(toExcl)), trind, trelem)));
                        dat.tri = arrayfun(@(x,y) dat.tri(x,y)-sum(toExcl(1:dat.tri(x,y))), trind, trelem)';
                        dat.tri(indtri,:) = [];
                    end
                    allPos(:,:,i) = dat.pos;
                    allInp{i} = dat;
                end
                for i = 1:numel(allInp)
                    allInp{i}.pos = mean(allPos,3);
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
            if ~isfield(allInp{1},'time')
                cfg{1}.latency = 'all';
                cfg{1}.correctm = thr.correctiontimepoint;
            end
            switch inpType
                case 'timefreq' % separate analyses for each band
                    fwoi = aas_getsetting(aap,'diagnostics.snapshotfwoi');
                    if ~isempty(fwoi)
                        for f = 1:size(fwoi,1)
                            cfg{f} = cfg{1};
                            cfg{f}.frequency = fwoi(f,:);
                            savepath{f} = savepath{1};
                        end
                    else
                        cfg{1}.frequency = 'all';
                    end                    
                case 'peak' % custom analysis and plotting
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
                avgcfg = keepfields(cfg{c},{'parameter','latency'});
                if isfield(pcfg,'snapshottwoi') && ~isempty(pcfg.snapshottwoi)
                    pcfg.snapshottwoi = pcfg.snapshottwoi(...
                        pcfg.snapshottwoi(:,1)/1000 >= stat.time(1) & ...
                        pcfg.snapshottwoi(:,2)/1000 <= stat.time(end) ...
                        ,:);
                end
                
                if numel(unique(m.groupmodel)) <= 2
                    for g = unique(m.groupmodel)
                        groupStat{end+1} = ft_granddescriptives(avgcfg, allInp{cfg{c}.design(1,:)==g});
                    end
                else
                    groupStat{1} = ft_granddescriptives(avgcfg, allInp{:});
                end
                if isfield(cfg{c},'frequency')
                    groupStat = cellfun(@(x) ft_selectdata(struct('frequency',cfg{c}.frequency,'avgoverfreq','yes'),x), groupStat,'UniformOutput',false);
                    pcfg.snapshotfwoi = cfg{c}.frequency;
                    cfg{c}.correctm = sprintf('%s_freq-%1.2f-%1.2f',cfg{c}.correctm,cfg{c}.frequency); % hack to add suffix to statFn
                end
                groupStat{1}.stat = stat;
                pcfg.parameter = cfg{c}.parameter;
                fdiag(groupStat,pcfg,m.name,statFn);
                            
                groupStat{1}.stat.subjects = m.subjects(subjmodel);
                statFn = fullfile(aas_getstudypath(aap),[m.name '_' cfg{c}.correctm '.mat']);
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