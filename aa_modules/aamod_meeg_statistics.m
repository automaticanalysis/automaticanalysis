
function [aap, resp] = aamod_meeg_statistics(aap,task)

resp='';

SETTING2DIM = containers.Map(...
    {'snapshotfwoiphase' 'snapshotfwoiamplitude' 'snapshotphwoi' 'snapshottwoi'},...
    {'freqlow' 'freqhigh' 'phase' 'time'}...
    );
DIM2SETTING = containers.Map(...
    {'freqlow' 'freqhigh' 'phase' 'time'},...
    {'snapshotfwoiphase' 'snapshotfwoiamplitude' 'snapshotphwoi' 'snapshottwoi'}...
    );

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
        [~, FT] = aas_cache_get(aap,'fieldtrip');
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
        if strcmp(inpstreams{1},'crossfreq') % tweak data into TFR
            data.label = data.elec.label; 
            data.freq = data.freqlow;
            data.time = data.freqhigh;
            data.dimord = 'chancmb_freq_time';
        end
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
            case 'crossfreq'
                statcfg.parameter   = 'crsspctrm';
                fstat = @meeg_statistics;
                fdiag = @meeg_diagnostics_CF;
                thr.correctiontimeseries = thr.correction;
                thr.correctiontimepoint = thr.correction;
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
        
        hashighressurface = all(arrayfun(@(subj) aas_stream_has_contents(aap,'subject',subj,'sourcesurface'), 1:aas_getN_bydomain(aap,'subject')));
        if hashighressurface && ~strcmp(statplotcfg.surface,'sourcemodel')
            fnsurf = cellstr(aas_getfiles_bystream(aap,'subject',1,'sourcesurface'));
            fnsurf = fnsurf{contains(fnsurf,statplotcfg.surface)};            
            grouphighressurface = ft_read_headshape(fnsurf);
            for subj = 2:aas_getN_bydomain(aap,'subject')
                fnsurf = cellstr(aas_getfiles_bystream(aap,'subject',subj,'sourcesurface'));
                fnsurf = fnsurf{contains(fnsurf,statplotcfg.surface)};
                mesh = ft_read_headshape(fnsurf);
                grouphighressurface.pos = grouphighressurface.pos + mesh.pos;
            end
            grouphighressurface.pos = grouphighressurface.pos/aas_getN_bydomain(aap,'subject');
            statplotcfg.surface = grouphighressurface;
        else
            if isfield(statplotcfg,'surface'), statplotcfg = rmfield(statplotcfg,'surface'); end
        end
        
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
                    aas_log(aap,false,sprintf('INFO: %d samples selected\n\tThe final subject list is saved in groupStat{1}.stat.subjects.',i));
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
                    [dat,perm,nshift] = shiftdata(dat,dimTime);
                    dat = dat(allTime,:,:,:,:);
                    dat = squeeze(unshiftdata(dat,perm,nshift));
                    allInp{i} = setfield(allInp{i}, par{:}, dat);
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
            
            if isfield(allInp{1},'labelcmb')
                aas_log(aap,false,'INFO: selecting common channelcombinations');
                labelcmb = categorical(allInp{1}.labelcmb);
                for subj = 2:numel(allInp)
                    labelcmb = intersect(labelcmb,categorical(allInp{subj}.labelcmb),'rows');
                end
                for subj = 1:numel(allInp)
                    [~,~,ind] = intersect(labelcmb,categorical(allInp{subj}.labelcmb),'rows');
                    allInp{subj}.labelcmb = allInp{subj}.labelcmb(ind,:);
                    allInp{subj}.crsspctrm = allInp{subj}.crsspctrm(ind,:,:,:,:);
                end
                
                % update neighbours for channelcombination
                cmblabels = cellstr(spm_file(num2str([1:size(labelcmb,1)]'),'prefix','cmb'));
                for lc1 = 1:size(labelcmb,1)
                    nb = false(size(labelcmb,1),1);
                    for lc2 = 1:size(labelcmb,2)
                        nb = nb | any(labelcmb(lc1,lc2) == labelcmb,2);
                        if thr.neighbours > 0, nb = nb |...
                                any(cell2mat(cellfun(@(l) l == labelcmb, neighbours(labelcmb(lc1,lc2) == {neighbours.label}).neighblabel, 'UniformOutput', false)),2);
                        end
                    end
                    cmbneighbours(lc1).label = cmblabels{lc1};
                    cmbneighbours(lc1).neighblabel = cmblabels(nb);
                end
                cfg{1}.neighbours = cmbneighbours;
                if cfg{1}.minnbchan == 0, cfg{1}.minnbchan = 2; end
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
                case 'crossfreq'
                    diag = aas_getsetting(aap,'diagnostics');
                    meas = fieldnames(diag);
                    numtask = cellfun(@(f) size(diag.(f),1), meas);
                    meas = meas(numtask>0);
                    numtask = numtask(numtask>0);
                    cntr = ones(1,numel(meas));
                    for c = 1:prod(numtask)
                        cfg{c} = cfg{1};
                        for indm = 1:numel(meas)
                            cfg{c}.average.(SETTING2DIM(meas{indm})) = diag.(meas{indm})(cntr(indm),:);
                            if strcmp(meas{indm},'snapshottwoi'), cfg{c}.average.time = cfg{c}.average.time./1000; end % correct time
                            savepath{c} = savepath{1};
                        end
                        cntr(1) = cntr(1) + 1;
                        for cc = 1:numel(cntr)-1
                            if cntr(cc) > size(diag.(meas{cc}),1)
                                cntr(cc) = 1;
                                cntr(cc+1) = cntr(cc+1) + 1;
                            end
                        end
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
                groupStat = cellfun(@(x) struct_update(x,allInp{1},'Mode','extend'), groupStat,'UniformOutput',false);
                if isfield(cfg{c},'frequency')
                    groupStat = cellfun(@(x) ft_selectdata(struct('frequency',cfg{c}.frequency,'avgoverfreq','yes'),x), groupStat,'UniformOutput',false);
                    if ischar(cfg{c}.frequency)
                        switch cfg{c}.frequency
                            case 'all'
                                cfg{c}.frequency = [min(groupStat{1}.freq) max(groupStat{1}.freq)];
                        end
                    end
                    pcfg.snapshotfwoi = cfg{c}.frequency;
                    cfg{c}.correctm = sprintf('%s_freq-%1.2f-%1.2f',cfg{c}.correctm,cfg{c}.frequency); % add suffix to statFn
                end
                if isfield(cfg{c}, 'average')
                    avgcfg = keepfields(cfg{c},{'parameter','average'});
                    avgcfg.parameter = {avgcfg.parameter [avgcfg.parameter 'sem']};
                    groupStat = cellfun(@(x) meeg_average(avgcfg,x), groupStat,'UniformOutput',false);
                    for fave = fieldnames(cfg{c}.average)'
                        pcfg.(DIM2SETTING(fave{1})) = cfg{c}.average.(fave{1});
                        statFn = sprintf('%s_%s-%1.2f-%1.2f',statFn,fave{1},cfg{c}.average.(fave{1})); % add suffix to statFn
                    end
                end
                groupStat{1}.stat = stat;
                pcfg.parameter = cfg{c}.parameter;
                fdiag(groupStat,pcfg,m.name,statFn);
                            
                groupStat{1}.stat.subjects = m.subjects(subjmodel);
                statFn = fullfile(aas_getstudypath(aap),[m.name '_' cfg{c}.correctm '.mat']);
                save(statFn,'groupStat');
                
                % append to stream
                outputFn = {};
                outstreamFn = aas_getoutputstreamfilename(aap,'study',[],'groupstat');
                if exist(outstreamFn,'file')
                    outputFn = cellstr(aas_getfiles_bystream(aap,'study',[],'groupstat','output'));
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

function data = meeg_average(cfg,data)
for task = fieldnames(cfg.average)'
    if ischar(cfg.average.(task{1}))
        switch cfg.average.(task{1})
            case 'all'
                cfg.average.(task{1}) = [min(data.(task{1})) max(data.(task{1}))];
        end
    end
    
    origdimord = strsplit(data.dimord,'_');
    taskdim = find(strcmp(origdimord,task{1}));
    taskind = data.(task{1}) >= cfg.average.(task{1})(1) & data.(task{1})<= cfg.average.(task{1})(2);
    
    if ~iscell(cfg.parameter), cfg.parameter = {cfg.parameter}; end
    for par = cfg.parameter
        dat = data.(par{1});
        [dat,perm,nshift] = shiftdata(dat,taskdim);
        dat = mean(dat(taskind,:,:,:,:),1);
        data.(par{1}) = squeeze(unshiftdata(dat,perm,nshift));
    end
    
    data.(task{1}) = mean(cfg.average.(task{1}));
    data.dimord = strjoin(setdiff(origdimord,task),'_');
end
end

function stat = meeg_statistics(cfg,varargin)
% this is a generic function without any sanitycheck
cfg.channel = {cfg.neighbours.label};

% average if needed
if isfield(cfg,'average')
    varargin = cellfun(@(x) meeg_average(cfg,x), varargin,'UniformOutput',false);
end

cfg.dim = size(varargin{1}.(cfg.parameter));
cfg.dimord = varargin{1}.dimord;

dat = cell2mat(cellfun(@(x) x.(cfg.parameter)(:), varargin, 'UniformOutput', false));

statmethod = str2func(['ft_statistics_' cfg.method]);
stat = statmethod(cfg, dat, cfg.design);

for fn = fieldnames(stat)'
    if size(stat.(fn{1}),1) == prod(cfg.dim)
        stat.(fn{1}) = reshape(stat.(fn{1}), cfg.dim);
    end
end

stat.cfg = cfg;
stat = struct_update(stat,varargin{1},'Mode','extend');

end