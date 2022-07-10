function [aap, resp] = aamod_meeg_statistics(aap,task)

resp='';

SETTING2DIM = containers.Map(...
    {'snapshotfwoiphase' 'snapshotfwoiamplitude' 'snapshotphwoi' 'snapshotfwoi'},...
    {'freqlow' 'freqhigh' 'phase' 'freq'}...
    );
DIM2SETTING = containers.Map(...
    {'freqlow' 'freqhigh' 'phase' 'freq'},...
    {'snapshotfwoiphase' 'snapshotfwoiamplitude' 'snapshotphwoi' 'snapshotfwoi'}...
    );

SETTING2CFG = containers.Map(...
    {'snapshotfwoi'},...
    {'frequency'}...
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
        
        aap = aas_report_add(aap,[],'<table id="data"><tr>');
        for m = models, aap = aas_report_add(aap,[],['<th>Model: ' m.name '</th>']); end
        aap = aas_report_add(aap,[],'</tr><tr>');
        
        for m = models
            aap = aas_report_add(aap,[],'<td valign="top">');
            
            savepath = fullfile(aas_getstudypath(aap),['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_' m.name]);            
            res = cellstr(spm_select('FPList',spm_file(savepath,'path'),['^' spm_file(savepath,'basename') '_.*']));
            if isempty(res{1}), continue; end
            % sort images according to datetime
            resd = cellfun(@dir, res);
            [~,ind] = sort(datenum({resd.date}));
            res = res(ind);
            
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
        dat = load(allFnTL{1}{1}); data = dat.(char(fieldnames(dat))); 
        cfg = [];        
        if isfield(data,'elec')
            cfg.elec = data.elec;
            
            if ~isempty(aas_getsetting(aap,'neighbourhood')) && strcmp(aas_getsetting(aap,'neighbourhood.method'),'distance')
                cfg.feedback    = 'no';
                cfg.method = 'distance';
                if ischar(aas_getsetting(aap,'neighbourhood.distance'))
                    switch aas_getsetting(aap,'neighbourhood.distance')
                        case 'autosymm'
                            % estimate distance between 1-100 mm
                            % to achieve symmetric neighbourhood
                            % while favouring more neighbours
                            cost = Inf(1,100);
                            for dist = 1:100
                                cfg.neighbourdist = dist;
                                neigh = ft_prepare_neighbours(cfg);
                                if all(arrayfun(@(n) numel(n.neighblabel), neigh) > 0)
                                    cost(dist) = sum(arrayfun(@(r) abs(numel(neigh(r).neighblabel)-numel(neigh(r+31).neighblabel)), 1:31))/(mean(arrayfun(@(n) numel(n.neighblabel), neigh))^(1/3));
                                end
                            end
                            [~,dist] = sort(cost);
                            cfg.neighbourdist = dist(1);
                    otherwise % method:parameter
                        method = strsplit(aas_getsetting(aap,'neighbourhood.distance'),':');
                        [method, parameter] = deal(method{:}); parameter = str2double(parameter);
                        switch method
                            case 'minimumneighbour'
                                for dist = 1:100
                                    cfg.neighbourdist = dist;
                                    neigh = ft_prepare_neighbours(cfg);
                                    if all(arrayfun(@(n) numel(n.neighblabel), neigh) >= parameter), break; end
                                end
                                cfg.neighbourdist = dist;
                        end
                    end
                else
                    cfg.neighbourdist = aas_getsetting(aap,'neighbourhood.distance');
                end
            else
                cfg.method      = 'triangulation';
            end
            
            f = figure;
            if exist('cost','var')
                set(f,'position',[0,0,1440 720]);
                subplot(1,2,1);
                plot(cost);
                p = subplot(1,2,2);
            else
                set(f,'position',[0,0,720 720]);
                p = subplot(1,1,1);
            end
            set(f,'PaperPositionMode','auto');
            cfg.feedback    = 'yes';
            neighbours = ft_prepare_neighbours(cfg);
            currfig = gcf;
            arrayfun(@(x) copyobj(x,p), flipud(get(get(currfig,'CurrentAxes'),'Children')));
            close(currfig);
            axis equal            
            
            print(f,'-noui',fullfile(aas_getstudypath(aap),['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_neighbours']),'-djpeg','-r300');
            close(f);
        else
            neighbours = [];
        end
        
        statplotcfg = aas_getsetting(aap,'diagnostics');
        if isfield(data,'elec')
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
            case {'timefreq' 'timeband'}
                statcfg.parameter   = 'powspctrm';
                fstat = @ft_freqstatistics;
                fdiag = @meeg_diagnostics_TFR;
            case {'crossfreq' 'crossband'}
                statcfg.parameter   = 'crsspctrm';
                fstat = @meeg_statistics;
                fdiag = @meeg_diagnostics_conn;
            case {'connfreq' 'connband'}
                f = fieldnames(data);
                f = f(cellfun(@(x) ~isempty(regexp(x,'.*spctrm$', 'once')),f));
                if numel(f) ~= 1, aas_log(aap,true,'Parameter cannot be identified'); end
                statcfg.parameter   = f{1};
                fstat = @meeg_statistics;
                fdiag = @meeg_diagnostics_conn;
        end
        if ft_datatype(data,'source')
            switch inpstreams{1}
            case {'timefreq' 'timeband'}
                statcfg.parameter   = 'pow';
            end
            fstat = @ft_sourcestatistics;
            fdiag = @meeg_diagnostics_source;
        end
        statcfg.alpha       = thr.p;
        statcfg.tail        = 0; % two-tailed
        statcfg.correcttail = 'prob';
        statcfg.method              = thr.method;
        statcfg.correctm            = thr.correction;
        statcfg.clusteralpha        = thr.p;
        statcfg.numrandomization    = thr.iteration;
        statcfg.minnbchan           = thr.neighbours;      % minimal number of neighbouring channels
        statcfg.neighbours          = neighbours; % defined as above
        statcfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable        
        
        if isfield(data,'tri') && ~strcmp(statplotcfg.background,'sourcemodel') && ... % only for cortical sheet source-level data
                all(arrayfun(@(subj) aas_stream_has_contents(aap,'subject',subj,'sourcesurface'), 1:aas_getN_bydomain(aap,'subject'))) && ... % only when all subject has surface
                ~aas_stream_has_contents(aap,'study',[],'background') % only when it has not been generated, yet
            fnsurf = cellstr(aas_getfiles_bystream(aap,'subject',1,'sourcesurface'));
            fnsurf = fnsurf{contains(fnsurf,statplotcfg.background)};            
            grouphighressurface = ft_read_headshape(fnsurf);
            for subj = 2:aas_getN_bydomain(aap,'subject')
                fnsurf = cellstr(aas_getfiles_bystream(aap,'subject',subj,'sourcesurface'));
                fnsurf = fnsurf{contains(fnsurf,statplotcfg.background)};
                mesh = ft_read_headshape(fnsurf);
                grouphighressurface.pos = grouphighressurface.pos + mesh.pos;
            end
            grouphighressurface.pos = grouphighressurface.pos/aas_getN_bydomain(aap,'subject');
            statplotcfg.background = grouphighressurface;
            outputFn = fullfile(aas_getstudypath(aap),'background.mat');
            save(outputFn,'grouphighressurface');
            aap = aas_desc_outputs(aap,'study',[],'background',outputFn);
        elseif isfield(data,'dim')
            switch statplotcfg.background
                case 'template'
                    dat = load(fullfile(FT.toolPath, 'template/headmodel/standard_mri.mat'));
                    statplotcfg.background = dat.mri;
                otherwise % filename
                    statplotcfg.background = ft_read_mri(statplotcfg.background);
            end
        else
            if isfield(statplotcfg,'background'), statplotcfg = rmfield(statplotcfg,'background'); end
        end
                
        atlascfg = [];
        atlascfg.interpmethod = 'nearest';
        atlascfg.parameter = 'aparc';
        
        diag = aas_getsetting(aap,'diagnostics');
        fieldwoi = {}; fieldofinterest = {}; indDiag = {}; indData = {};
        if isfield(diag,'snapshotfwoi'), fieldwoi{end+1} = 'snapshotfwoi'; fieldofinterest{end+1} = 'band'; end
        if isfield(diag,'snapshotfwoiphase') && ~isempty(diag.snapshotfwoiphase), fieldwoi{end+1} = 'snapshotfwoiphase'; fieldofinterest{end+1} = 'bandlow'; end
        if isfield(diag,'snapshotfwoiamplitude') && ~isempty(diag.snapshotfwoiamplitude), fieldwoi{end+1} = 'snapshotfwoiamplitude'; fieldofinterest{end+1} = 'bandhigh'; end
        if isempty(fieldwoi), aas_log([],true,'no valid snapshot specification found'); end
        
        models = aas_getsetting(aap,'model'); models(1) = [];
        % do not redo fully analyzed models
        nRes = prod(cellfun(@(f) size(diag.(f),1), fieldwoi));
        modelIsRun = arrayfun(@(m) ...
            numel(cellstr(spm_select('List',aas_getstudypath(aap),['^' m.name '_']))) == nRes |...
            numel(cellstr(spm_select('List',aas_getstudypath(aap),['^' m.name]))) == nRes,...
            models);
        if ~isempty(modelIsRun)
            models(modelIsRun) = []; 
            % clean existig results for unfinished models
            if aas_stream_has_contents(aap,'study',[],'groupstat')
                resStat = cellstr(aas_getfiles_bystream(aap,'study',[],'groupstat'));
                resFns = arrayfun(@(m) cellstr(spm_select('FPList',aas_getstudypath(aap),['^' m.name])),models,'UniformOutput',false);
                for m = cat(1,resFns{:})'
                    delete(m{1});
                    resStat(strcmp(resStat,m{1})) = [];
                end
                aap = aas_desc_outputs(aap,'study',[],'groupstat',resStat);
            end
        end
        outputFn = {};
        for m = models
            diag = aas_getsetting(aap,'diagnostics');
            atlas = [];
            savepath = fullfile(aas_getstudypath(aap),['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_' m.name]);
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
              
            %% Retrieve Data
            allInp = {}; subjmodel = [];
            for subj = 1:numel(m.subjects)
                aasubj = find(strcmp({aap.acq_details.subjects.subjname},m.subjects{subj}));
                if ~any(m.trialmodel{1}=='_') % ER/TFR
                    inpType = inpstreams{1};
                    sFn = allFnTL{aasubj};
                else % peak
                    inpType = 'peak';
                    sFn = allFnP{aasubj};
                end                
                for iFn = find(sum(cell2mat(cellfun(@(y) cellfun(@(x) ~isempty(x), regexp(spm_file(sFn,'basename'),[inpType '_' y '$'])), m.trialmodel,'UniformOutput',false)),2))
                    if isempty(iFn), aas_log(aap,true,'Trialmodels not found'); end
                    for fn = sFn(iFn)'
                        aas_log(aap,false,sprintf('INFO: load %s',fn{1}));
                        dat = load(fn{1}); dat = dat.(string(fieldnames(dat)));
                        
                        % band -> freq/diag.snapshotfwoi
                        if any(isfield(dat,{'band','bandlow','bandhigh'})) % assume same data type
                            for f = 1:numel(fieldwoi)
                                [~,indDiag{f},indData{f}] = intersect(diag.(fieldwoi{f}),dat.(fieldofinterest{f}),'stable');
                                dat = rmfield(dat,intersect(fieldnames(dat),{fieldofinterest{f} 'bandspec'}));
                                dat.(strrep(fieldofinterest{f},'band','freq')) = reshape(indData{f},1,[]);
                                if isfield(dat,'avg')
                                    indDim = find(strcmp(strsplit(dat.avg.dimord,'_'),fieldofinterest{f}));
                                    [meas,perm,nshift] = shiftdata(dat.avg.(statcfg.parameter),indDim);
                                    meas = meas(indData{f},:,:,:,:);
                                    dat.avg.(statcfg.parameter) = unshiftdata(meas,perm,nshift);
                                    dat.avg.dimord = strrep(dat.avg.dimord,'band','freq');
                                else
                                    indDim = find(strcmp(strsplit(dat.dimord,'_'),fieldofinterest{f}));
                                    [meas,perm,nshift] = shiftdata(dat.(statcfg.parameter),indDim);
                                    meas = meas(indData{f},:,:,:,:);
                                    dat.(statcfg.parameter) = unshiftdata(meas,perm,nshift);
                                    dat.dimord = strrep(dat.dimord,fieldofinterest{f},strrep(fieldofinterest{f},'band','freq'));
                                end
                            end
                        end

                        if strcmp(cfg{1}.channel, 'aa') % split data into channels to be able to model them
                            chdata = {};
                            for ch = m.channels.channels
                                chdata{end+1} = ft_selectdata(struct('channel',ch),dat);
                                chdata{end}.label = cfg{1}.channel;
                                chdata{end}.elec.label(strcmp(chdata{end}.elec.label,ch)) = cfg{1}.channel;
                            end
                            dat = ft_combine(combinecfg,chdata{:});
                        end                        
                        % expand data and remove confouding fields
                        if ~isfield(dat,statcfg.parameter)
                            dat = ft_selectdata([],dat);
                            dat = rmfield(dat,intersect(fieldnames(dat),{'ori','filterdimord','inside'}));
                            if all(isfield(dat,{'pos' 'label'})), dat = rmfield(dat,'label'); end
                            dat.cfg = keepfields(dat.cfg.previous,'included');
                        end
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
                        if isfield(dat,'dim') % grid-based source -> replace and extend pos with that of the template
                            tmpdat = load(aas_getfiles_bystream(aas_setcurrenttask(aap,aas_getsourcestage(aap,'aamod_meeg_sourcereconstruction')),'subject',subj,'sourcemodel'));
                            tmpdat = load(tmpdat.sourcemodel.cfg.template);
                            dat.pos = tmpdat.sourcemodel.pos;
                            parsize = size(dat.(statcfg.parameter)); parsize(1) = numel(dat.cfg.included);
                            par = zeros(parsize); par(dat.cfg.included,:,:,:,:) = dat.(statcfg.parameter);
                            dat.(statcfg.parameter) = par;
                        end
                        dat.cfg = keepfields(dat.cfg,'included'); % remove provenance to save space
                        allInp{end+1} = ft_struct2single(dat);
                        subjmodel(end+1) = subj;
                    end
                end
            end
            if exist('indDiag','var')
                diagBands = {};
                for f = 1:numel(fieldwoi)
                    diagBands{f} = diag.(fieldwoi{f});
                    diag.(fieldwoi{f}) = [indDiag{f}-0.25 indDiag{f}+0.25];
                    pcfg.(fieldwoi{f}) = [indDiag{f}-0.25 indDiag{f}+0.25];
                end
            end
            
            %% Check temporal overlap (for time-resolved data)
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
                            i = find(atrial==max(atrial),1,'first');
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
                if ~isfield(allInp{1},'dimord'), aas_log(aap,true,'no dimension information found'); end
                dimord = allInp{1}.dimord;
                par = {statcfg.parameter};
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
                
                if strcmp(aas_getsetting(aap,'selectoverlappingdata.time'),'ignore') && ~isempty(pcfg.snapshottwoi) % equally divide trials without overlap
                    step = (max(allInp{i}.time) - min(allInp{i}.time))/(size(pcfg.snapshottwoi,1)-1);
                    pcfg.snapshottwoi = ((min(allInp{i}.time):step:max(allInp{i}.time))'-step/4)*1000;
                    pcfg.snapshottwoi(:,2) = pcfg.snapshottwoi(:,1)+step/2*1000;
                end
            end
            
            %% Check spatial overlap (for source)
            if isfield(allInp{1},'pos') % source
                allAtlas = struct('aparc',{},'aparclabel',{}, 'aparcpos', {});
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
                    if isfield(dat,'tri') % cortical sheet
                        dat.pos(toExcl,:) = [];
                        dat.(cfg{1}.parameter)(toExcl,:,:) = [];
                        [trind, trelem] = ndgrid(1:size(dat.tri,1),1:size(dat.tri,2)); trind = transpose(trind); trelem = transpose(trelem);
                        indtri = logical(sum(arrayfun(@(x,y) any(dat.tri(x,y)==find(toExcl)), trind, trelem)));
                        dat.tri = arrayfun(@(x,y) dat.tri(x,y)-sum(toExcl(1:dat.tri(x,y))), trind, trelem)';
                        dat.tri(indtri,:) = [];
                        allPos(:,:,i) = dat.pos;
                    end
                    allInp{i} = dat;
                    if aas_stream_has_contents(aap,'subject',subjmodel(i),'sourceatlas')
                        dat = load(aas_getfiles_bystream(aap,'subject',subjmodel(i),'sourceatlas')); sourceatlas = dat.sourceatlas;
                        dat = keepfields(ft_sourceinterpolate(atlascfg, sourceatlas, allInp{i}),fieldnames(allAtlas));
                        dat.aparcpos = sourceatlas.aparcpos;
                        allAtlas(i) = dat;
                    end
                end
                if isfield(allInp{1},'tri')
                    % make sure that all nodes are included in tri
                    allNodes = false(size(allInp{1}.pos,1),1);
                    tri = allInp{1}.tri; i = 1;
                    allNodes(unique(tri)) = true; i = 1;
                    while ~all(allNodes)
                        i = i + 1;
                        for n = intersect(allInp{i}.tri(:),find(~allNodes))'
                            tri = [tri; ...
                                allInp{i}.tri(find(arrayfun(@(t) any(allInp{i}.tri(t,:)==n),1:size(allInp{i}.tri,1)),1,'first'),:)...
                                ];
                            allNodes(unique(tri)) = true;
                        end
                    end
                    for i = 1:numel(allInp)
                        allInp{i}.pos = mean(allPos,3);
                        allInp{i}.tri = tri;
                    end
                    cfg{1}.tri = tri;
                end
                if ~isempty(allAtlas)
                    % check consitency
                    % - same list of areas
                    allParcLabMatch = sum(reshape(...
                        arrayfun(@(i1,i2) isequal(allAtlas(i1).aparclabel,allAtlas(i2).aparclabel),...
                            reshape(repmat(1:numel(allAtlas),numel(allAtlas),1),1,[]), reshape(repmat(1:numel(allAtlas),numel(allAtlas),1)',1,[])...
                            ),...
                        numel(allAtlas),numel(allAtlas)));
                    jointParcLab = allAtlas(find(allParcLabMatch==max(allParcLabMatch),1,'first')).aparclabel;
                    % - update parcellation
                    for i = 1:numel(allInp)
                        if ~isequal(allAtlas(i).aparclabel,jointParcLab)
                            % - extra area -> remove
                            missingParcLab = cellfun(@(a) ~any(strcmp(jointParcLab,a)), allAtlas(i).aparclabel);
                            adjParcLab = cumsum(missingParcLab); adjParcLab(missingParcLab) = -1;
                            if any(adjParcLab)
                                for a = 1:numel(adjParcLab)
                                    if adjParcLab(a) == -1, allAtlas(i).aparc(allAtlas(i).aparc == a) = 0;
                                    else
                                        allAtlas(i).aparc(allAtlas(i).aparc == a) = allAtlas(i).aparc(allAtlas(i).aparc == a) - adjParcLab(a);
                                    end
                                end
                                allAtlas(i).aparclabel = allAtlas(i).aparclabel(~missingParcLab);
                            end
                            
                            matchParcLab = cellfun(@(a) find(strcmp(jointParcLab,a)), allAtlas(i).aparclabel);
                            for a = 1:numel(matchParcLab)
                                if matchParcLab(a) ~= a
                                    allAtlas(i).aparc(allAtlas(i).aparc == a) = matchParcLab(a);
                                    allAtlas(i).aparclabel(a) = jointParcLab(a);
                                end
                            end
                        end
                    end
                    
                    % select dominant parcellation
                    atlas.aparc = mode([allAtlas.aparc],2);
                    atlas.aparclabel = jointParcLab;
                end                
            end
            
            %% Cluster-forming rule for channel combinations (connectivity and crossfrequency)
            if isfield(allInp{1},'labelcmb') && ~isfield(allInp{1},'label')
                aas_log(aap,false,'INFO: selecting common channelcombinations');
                labelcmb = categorical(allInp{1}.labelcmb);
                for subj = 2:numel(allInp)
                    labelcmb = intersect(labelcmb,categorical(allInp{subj}.labelcmb),'rows');
                end
                for subj = 1:numel(allInp)
                    [~,~,ind] = intersect(labelcmb,categorical(allInp{subj}.labelcmb),'rows');
                    allInp{subj}.labelcmb = allInp{subj}.labelcmb(ind,:);
                    allInp{subj}.(statcfg.parameter) = allInp{subj}.(statcfg.parameter)(ind,:,:,:,:);
                end
                
                % update neighbours for channelcombination
                cmblabels = cellstr(spm_file(num2str([1:size(labelcmb,1)]'),'prefix','cmb'));
                cmbneighbours = [];
                for lc1 = 1:size(labelcmb,1)
                    nb = false(size(labelcmb,1),1);
                    if ~thr.combinationneighbours
                        % option 1: do not consider channelcombination
                        if thr.neighbours > 0, nb = nb |...
                                (any(cell2mat(cellfun(@(l) labelcmb(:,2) == l, neighbours({neighbours.label} == labelcmb(lc1,2)).neighblabel, 'UniformOutput', false)),2) & labelcmb(:,1) == labelcmb(lc1,1)) |... % source to destinationneighbours
                                (any(cell2mat(cellfun(@(l) labelcmb(:,1) == l, neighbours({neighbours.label} == labelcmb(lc1,1)).neighblabel, 'UniformOutput', false)),2) & labelcmb(:,2) == labelcmb(lc1,2)); % destination to sourceneighbours
                        end
                    else
                        % option 2: consider channelcombination
                        for lc2 = 1:size(labelcmb,2)
                            nb = nb | any(labelcmb(lc1,lc2) == labelcmb,2);
                            if thr.neighbours > 0, nb = nb |...
                                    any(cell2mat(cellfun(@(l) l == labelcmb, neighbours(labelcmb(lc1,lc2) == {neighbours.label}).neighblabel, 'UniformOutput', false)),2);
                            end
                        end
                    end
                    cmbneighbours(lc1).label = cmblabels{lc1};
                    cmbneighbours(lc1).neighblabel = cmblabels(nb);
                end
                cfg{1}.neighbours = cmbneighbours;
                if cfg{1}.minnbchan == 0, cfg{1}.minnbchan = 2; end
            end
            
            %% Design
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
            
            %% Stat cfg for averages
            if ~isfield(allInp{1},'time'), cfg{1}.latency = 'all'; end
            switch inpType
                case {'timefreq' 'timeband' 'connfreq' 'connband' 'crossfreq' 'crossband'}                    
                    meas = intersect(fieldnames(diag),SETTING2DIM.keys);
                    numtask = cellfun(@(f) size(diag.(f),1), meas);
                    meas = meas(numtask>0);
                    numtask = numtask(numtask>0);                    
                    cntr = ones(1,numel(meas));
                    if ~isempty(cntr)
                        for c = 1:prod(numtask)
                            cfg{c} = cfg{1};
                            for indm = 1:numel(meas)
                                cfg{c}.average.(SETTING2DIM(meas{indm})) = diag.(meas{indm})(cntr(indm),:);
                                if any(strcmp({'timefreq' 'timeband'},inpType))
                                    cfg{c}.(SETTING2CFG(meas{indm})) = diag.(meas{indm})(cntr(indm),:);
                                    cfg{c}.savepath = [savepath '_topoplot'];
                                else
                                    cfg{c}.savepath = savepath;
                                end
                            end
                            cntr(1) = cntr(1) + 1;
                            for cc = 1:numel(cntr)-1
                                if cntr(cc) > size(diag.(meas{cc}),1)
                                    cntr(cc) = 1;
                                    cntr(cc+1) = cntr(cc+1) + 1;
                                end
                            end
                        end
                    end
                case 'peak' % custom analysis and plotting
                    cfg{1}.parameter = 'amp';
                    cfg{1}.savepath = [savepath '_amp'];
                    cfg{2} = cfg{1};
                    cfg{2}.parameter = 'lat';
                    cfg{2}.savepath = [savepath '_lat'];
            end
            
            %% Run
            if ~exist('diagBands','var'), diagBands = {}; fieldofinterest = ''; end
            doParallel = aas_getsetting(aap,'numberofworkers') > 1;
            if doParallel
                aapoolprofile = strsplit(aap.directory_conventions.poolprofile,':'); poolprofile = aapoolprofile{1};
                if ~strcmp(aap.options.wheretoprocess,'qsub'), aas_log(aap,false,sprintf('WARNING: pool profile %s is not used via DCS/MPaS; therefore it may not work for parfor',poolprofile)); end
                try
                    cluster = parcluster(poolprofile);
                    if numel(aapoolprofile) > 1, cluster.SubmitArguments = aapoolprofile{2}; end
                    cluster.SubmitArguments = compose("%s --mem=%dG --time=%d", string(cluster.SubmitArguments), aap.options.aaparallel.memory, aap.options.aaparallel.walltime*60);
                    global aaworker;
                    wdir = spm_file(tempname,'basename'); wdir = fullfile(aaworker.parmpath,wdir(1:8));
                    aas_log(aap,false,['INFO: parallel job storage is ' wdir]);
                    aas_makedir(aap,wdir);
                    cluster.JobStorageLocation = wdir;
                    nWorkers = min([cluster.NumWorkers aas_getsetting(aap,'numberofworkers')]);
                catch E
                    aas_log(aap,false,['WARNING: ' poolprofile ' could not been initialised - ' E.message ' --> parallelisation is disabled']);
                    doParallel = false;
                end
            end
            if doParallel
                global aaworker
                maxIterPerJob = ceil(numel(cfg)/nWorkers);
                cfg_perm = cfg(1); cfg_perm(1) = []; % create empty array of cfgs
                for indCfg = 1:numel(cfg)
                    cfg_perm(end+1) = cfg(indCfg);
                    if (numel(cfg_perm) == maxIterPerJob) ||... % batch is ready
                            (indCfg == numel(cfg)) % last batch
                        batch(cluster,@run_model,1,{aap,m,fstat,fdiag,cfg_perm,pcfg,allInp,subjmodel,atlas,DIM2SETTING,fieldofinterest,diagBands,aaworker},'AutoAttachFiles',false); % run permutations
                        cfg_perm(1:end) = []; % reset array of cfgs
                    end
                end
                
                while ~all(arrayfun(@(j) strcmp(j.State,'finished'), cluster.Jobs)), pause(1); end
                statFn = arrayfun(@fetchOutputs, cluster.Jobs);
                statFn = cat(1,statFn{:});
            else
                statFn = run_model(aap,m,fstat,fdiag,cfg,pcfg,allInp,subjmodel,atlas,DIM2SETTING,fieldofinterest,diagBands,[]);
            end
            
            % collect stasts
            outputFn = cat(1,outputFn,cellstr(statFn));
        end
        
        % append to stream
        outstreamFn = aas_getoutputstreamfilename(aap,'study',[],'groupstat');
        if exist(outstreamFn,'file')
            outputFn0 = cellstr(aas_getfiles_bystream(aap,'study',[],'groupstat','output'));
        end
        outputFn = cat(1,outputFn0,outputFn);
        aap = aas_desc_outputs(aap,'study',[],'groupstat',outputFn);

        FT.unload;

    case 'checkrequirements'
        if ~aas_cache_get(aap,'fieldtrip'), aas_log(aap,true,'FieldTrip is not found'); end
        
        % idenitfy use case and adjust other inputs
        inpstreams = aas_getstreams(aap,'input');
        % - not source
        tmpaap = aas_findinputstreamsources(aap,false); % update streams
        tmpaap.options.verbose = -1;
        if isempty(aas_getsourcestage(tmpaap,'aamod_meeg_sourcereconstruction',inpstreams{1}))
            if aas_stream_has_contents(aap,'sourcesurface'), aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'sourcesurface',[]); end
            if aas_stream_has_contents(aap,'sourceatlas'), aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'sourceatlas',[]); end
        end

        % - band
        if contains(inpstreams{1},'band')
            diag = aas_getsetting(aap,'diagnostics'); dg = fieldnames(diag);
            for fw = dg(contains(dg,'fwoi'))'
                if ~isempty(diag.(fw{1})) && ~ischar(diag.(fw{1}){1})
                    aas_log(aap,true,'For band-type input, diagnostics.snapshotfwoi* MUST be a subset of the bandnames');
                end
            end
        end
end
end

function statFn = run_model(aap,m,fstat,fdiag,cfg,pcfg,allInp,subjmodel,atlas,DIM2SETTING,fieldofinterest,diagBands,aaworker)
if isstruct(aaworker) && isfield(aaworker,'aacache')
    global aacache
    aacache = aaworker.aacache;
end

for c = 1:numel(cfg)
    diagFn = cfg{c}.savepath;
    
    stat = fstat(cfg{c}, allInp{:});
    
    % plot
    groupStat = {};
    avgcfg = keepfields(cfg{c},{'parameter','latency'});
    if isfield(pcfg,'snapshottwoi') && ~isempty(pcfg.snapshottwoi)
        pcfg.snapshottwoi = pcfg.snapshottwoi(...
            arrayfun(@(i) any(arrayfun(@(t) pcfg.snapshottwoi(i,1) < t & pcfg.snapshottwoi(i,2) > t, stat.time*1000)), 1:size(pcfg.snapshottwoi,1)) ...
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
    if isfield(cfg{c}, 'average')
        avgcfg = keepfields(cfg{c},{'parameter','average'});
        avgcfg.parameter = {avgcfg.parameter [avgcfg.parameter 'sem']};
        groupStat = cellfun(@(x) meeg_average(avgcfg,x), groupStat,'UniformOutput',false);
        for fave = fieldnames(cfg{c}.average)'
            pcfg.(DIM2SETTING(fave{1})) = cfg{c}.average.(fave{1});
            if contains(fave{1},'freq') && exist('diagBands','var')
                indMeas = strcmp(fieldofinterest,strrep(fave{1},'freq','band'));
                diagFn = sprintf('%s_%s-%s',diagFn,fave{1},strrep(diagBands{indMeas}{mean(cfg{c}.average.(fave{1}))},' ','')); % add suffix to statFn
            else
                diagFn = sprintf('%s_%s-%1.2f-%1.2f',diagFn,fave{1},cfg{c}.average.(fave{1})); % add suffix to statFn
            end
        end
    end
    stat.mask(isnan(stat.mask)) = 0;
    groupStat{1}.stat = stat;
    groupStat{1}.stat.subjects = m.subjects(subjmodel);
    if isstruct(atlas), groupStat{1} = struct_update(groupStat{1},atlas,'Mode','extend'); end
    statFn{c} = spm_file(strrep(diagFn,['diagnostic_' aap.tasklist.main.module(aap.tasklist.currenttask.modulenumber).name '_'],''),'ext','mat');
    save(statFn{c},'groupStat');
    
    pcfg.parameter = cfg{c}.parameter;
    fdiag(groupStat,pcfg,m.name,diagFn);
end
statFn = reshape(statFn,[],1);
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