function h = meeg_diagnostics_ER(er,diag,figtitle,savepath)
TOPO_MAXCOL = 4; % maximum number of columns of topoplot mosaic 
TOL = 1e-8;

if nargin < 3, figtitle = 'Sample'; end

if ~iscell(er), er = {er}; end

if ~isfield(diag,'layout'), diag.layout = ft_prepare_layout([],er{1}); end
dispcfg = [];
dispcfg.parameter = diag.parameter;
dispcfg.layout = diag.layout;
dispcfg.marker = 'labels';
dispcfg.interactive = 'no';
dispcfg.showoutline = 'yes';
dispcfg.comment = 'no';

cfgmulti = dispcfg;
cfgtopo = dispcfg;
cfgtopo.figure = 'gca';

cfgmulti.showlabels  = 'yes';
if isfield(er{1},'stat')
    stat = er{1}.stat;
    er{1} = rmfield(er{1},'stat');
    if strcmp(stat.cfg.correctm,'cluster')
        if isfield(stat,'posclusters')
            pos_cluster_pvals = [stat.posclusters(:).prob];
            pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
            pos = ismember(stat.posclusterslabelmat, pos_signif_clust);
        else
            pos = stat.mask*0;
        end
        
        if isfield(stat,'posclusters')
            neg_cluster_pvals = [stat.negclusters(:).prob];
            neg_signif_clust = find(neg_cluster_pvals < stat.cfg.alpha);
            neg = ismember(stat.negclusterslabelmat, neg_signif_clust);
        else
            neg = stat.mask*0;
        end
    else
        pos = stat.mask & (stat.stat > 0);
        neg = stat.mask & (stat.stat < 0);
    end
    er{1}.mask = stat.mask;
    cfgmulti.maskparameter = 'mask';
    cfgmulti.maskstyle = 'box';
    cfgtopo.highlight = 'on';
    [iLabER,iLabStat] = match_str(er{1}.label, stat.label);
    if isfield(diag,'topohighlight')
        topoHighlightFun = str2func(diag.topohighlight);
    else
        topoHighlightFun = @all;
    end
end

% add dummy time if none exist
if ~isfield(er{1},'time')
    for e = 1:numel(er)
        er{e}.time = 0;
    end
end

if numel(er{1}.time) > 1
    if numel(er{1}.label) == 1 % single channel data
        cfgmulti.linewidth = 2;
        ft_singleplotER(cfgmulti, er{:});
        diag.videotwoi =  [];
        diag.snapshottwoi = [];
    else
        ft_multiplotER(cfgmulti, er{:});
    end
    set(gcf,'Name',figtitle);
    if nargin >= 4 
        set(gcf,'Position',[0,0,1080 1080]);
        set(gcf,'PaperPositionMode','auto');
        print(gcf,'-noui',[savepath '_multiplot.jpg'],'-djpeg','-r300');
        close(gcf);
    else
        h.multiplot = gcf;
    end
else
    diag.videotwoi =  [];
    diag.snapshottwoi = [er{1}.time er{1}.time]*1000;
end

if numel(er) > 1
    labels = arrayfun(@(x) sprintf('group #%d',x),1:numel(er),'UniformOutput',false);
else
    labels = {''};
end

if exist('stat','var')
    er = [er(1) er];
    er{1}.(diag.parameter) = stat.stat;
    labels = [{'stat'} labels];
end
    
if ~isempty(diag.snapshottwoi)
    indER = 1:numel(er);
    if strcmp(labels{1},'stat'), indER = 2:numel(er); end
    if numel(er{1}.time) > 1
        zlimSmooth = mean(diff(diag.snapshottwoi'))/1000*1/mean(diff(er{1}.time));    
        alldat = cellfun(@(x) x.(diag.parameter), er(indER), 'UniformOutput', false);
        alldat = vertcat(alldat{:});
        alldat = arrayfun(@(x) smooth(alldat(x,:),zlimSmooth),1:size(alldat,1),'UniformOutput',false);
        alldat = horzcat(alldat{:});
        zlim = [min(alldat(:)) max(alldat(:))];
    else
        zlim = [min(cellfun(@(x) min(x.(diag.parameter)(:)),er(indER))) max(cellfun(@(x) max(x.(diag.parameter)(:)),er(indER)))];
    end
    
    f = figure('Name',figtitle);
    colPlot = min([TOPO_MAXCOL numel(er)*ceil(sqrt(size(diag.snapshottwoi,1)))]);
    colPlot = floor(colPlot/numel(er))*numel(er);
    rowPlot = ceil(size(diag.snapshottwoi,1)*numel(er)/colPlot);
    if rowPlot > colPlot % figure will likely not fit to screen
        warning('Topoplot figure is too large to fit to screen -> It will be hidden.')
        set(f,'Visible','off');
    end 
    set(f,'Position',[0,0,1080 ceil(1080/colPlot*rowPlot)]);
    
    for t = 1:size(diag.snapshottwoi,1)
        cfgtopo.xlim = diag.snapshottwoi(t,:)/1000;
        
        if exist('stat','var')
            % Get the index of each significant channel
            [junk, sampStart] = min(abs(er{1}.time - cfgtopo.xlim(1)));
            if er{1}.time(sampStart)+TOL < cfgtopo.xlim(1), sampStart=sampStart+1; end
            [junk, sampEnd] = min(abs(er{1}.time - cfgtopo.xlim(2)));
            if er{1}.time(sampEnd)-TOL > cfgtopo.xlim(2), sampEnd=sampEnd+1; end            
            
            pos_int = zeros(numel(er{1}.label),1);
            neg_int = zeros(numel(er{1}.label),1);
            pos_int(iLabER) = topoHighlightFun(pos(iLabStat, sampStart:sampEnd), 2);
            neg_int(iLabER) = topoHighlightFun(neg(iLabStat, sampStart:sampEnd), 2);
            cfgtopo.highlightchannel = find(pos_int | neg_int);
        end
        
        for e = 1:numel(er)
            cfgtopo.colorbar = 'no';
            if strcmp(labels{e},'stat')
                cfgtopo.zlim = [min(er{e}.(diag.parameter)(:)) max(er{e}.(diag.parameter)(:))]; 
                if t == 1, cfgtopo.colorbar = 'West'; end
            else
                cfgtopo.zlim = zlim;
                if t == 1 && (isempty(labels{e}) || strcmp(labels{e},'group #1')), cfgtopo.colorbar = 'West'; end
            end
            if (cfgtopo.zlim(1) < 0) && (cfgtopo.zlim(2) > 0)
                r = cfgtopo.zlim(2)/-cfgtopo.zlim(1);
                if r > 1, cmaps{t,e} = [cmapcold(round(64/r)); cmaphot(64)];
                else, cmaps{t,e} = [cmapcold(64); cmaphot(round(r*64))];
                end
            elseif cfgtopo.zlim(1) < 0, cmaps{t,e} = cmapcold(64);
            else
                cmaps{t,e} = cmaphot(64);
            end
            ax(t,e) = subplot(rowPlot,colPlot,(t-1)*numel(er)+e);
            title(sprintf('%s %03d-%03d ms',labels{e},diag.snapshottwoi(t,:)));
            ft_topoplotER(cfgtopo, er{e});
        end
    end
    % Set colormaps
    for t = 1:size(diag.snapshottwoi,1)
        for e = 1:numel(er)
            if ~isempty(cmaps{t,e}), colormap(ax(t,e),cmaps{t,e}); end
        end
    end
    % Put colorbars between the axes
    cbs = findall(f,'type','ColorBar');
    axs = findall(f,'type','Axes');
    dPos = cbs(1).Position.*[0 0 0.25 0.25];
    dPos(1) = (-(diff(arrayfun(@(x) x.Position(1), axs(1:2))))-axs(1).Position(3))/2;    
    for c = cbs'
        c.YAxisLocation = 'left';
        c.Location = 'manual';
        c.Position = c.Position-dPos;
    end
    if nargin >= 4 
        set(f,'PaperPositionMode','auto');
        print(f,'-noui',[savepath '_topoplot.jpg'],'-djpeg','-r150');
        close(f)
    else
        h.topoplot = f;
    end
end
    
if ~isempty(diag.videotwoi) && (nargin >= 4) 
    dt = diag.videotwoi/1000;

    indER = 1:numel(er);
    if strcmp(labels{1},'stat'), indER = 2:numel(er); end
    zlimSmooth = dt*1/mean(diff(er{1}.time));
    alldat = cellfun(@(x) x.(diag.parameter), er(indER), 'UniformOutput', false);
    alldat = vertcat(alldat{:});
    alldat = arrayfun(@(x) smooth(alldat(x,:),zlimSmooth),1:size(alldat,1),'UniformOutput',false);
    alldat = horzcat(alldat{:});
    zlim = [min(alldat(:)) max(alldat(:))];
    cfgtopo.zlim = zlim;
    
    f = figure;
    set(f,'Position',[0,0,720 ceil(720/numel(er))]);
    
    aviObject = VideoWriter([savepath '_topoplot.avi']);
    aviObject.FrameRate = 1/dt;
    open(aviObject);
    
    t1 = er{1}.time(1);
    while t1+dt <= er{1}.time(end)
        clf;
        cfgtopo.xlim = [t1 t1+dt];
                
        if exist('stat','var')
            % Get the index of each significant channel
            [junk, sampStart] = min(abs(er{1}.time - cfgtopo.xlim(1)));
            if er{1}.time(sampStart)+TOL < cfgtopo.xlim(1), sampStart=sampStart+1; end
            [junk, sampEnd] = min(abs(er{1}.time - cfgtopo.xlim(2)));
            if er{1}.time(sampEnd)-TOL > cfgtopo.xlim(2), sampEnd=sampEnd+1; end            
            
            pos_int = zeros(numel(er{1}.label),1);
            neg_int = zeros(numel(er{1}.label),1);
            pos_int(iLabER) = topoHighlightFun(pos(iLabStat, sampStart:sampEnd), 2);
            neg_int(iLabER) = topoHighlightFun(neg(iLabStat, sampStart:sampEnd), 2);
            cfgtopo.highlightchannel = find(pos_int | neg_int);
        end
        
        for e = 1:numel(er)
            if e==1 || (e == 2 && strcmp(labels{1},'stat')), cfgtopo.colorbar = 'West';
            else, cfgtopo.colorbar = 'no'; end
            subplot(1,numel(er),e);
            ft_topoplotER(cfgtopo, er{e});
        end
        % Put colorbars between the axes
        cbs = findall(f,'type','ColorBar');
        axs = findall(f,'type','Axes');
        if numel(axs) == 1
            dPos = cbs(1).Position.*[0.25 0 0.25 0.25];
        else
            dPos = cbs(1).Position.*[0 0 0.25 0.25];
            dPos(1) = (-(diff(arrayfun(@(x) x.Position(1), axs(1:2))))-axs(1).Position(3))/2;
        end
        for c = cbs'
            c.YAxisLocation = 'left';
            c.Location = 'manual';
            c.Position = c.Position-dPos;
        end
        
        writeVideo(aviObject,getframe(f));
        t1 = t1+dt;
    end 
    
    close(aviObject);
    close(f);
end
end

%% Colormap functions
function cmap = cmapbase
[s, FT] = aas_cache_get([],'fieldtrip');
if s
    FT.load;
    FT.addExternal('brewermap');
    cmap = flipud(brewermap(128,'RdBu'));
else
    cmap = flipud([create_grad([0.4 0 0.1],[1 0 0],32);create_grad([1 0 0],[1 1 1],32);create_grad([1 1 1],[0 0 1],32);create_grad([0 0 1],[0 0.2 0.4],32)]);
end
end

function cmap = cmaphot(n)
cmap = cmapbase;
cmap = cmap(65:(end-(64-n)),:);
end

function cmap = cmapcold(n)
cmap = cmapbase;
cmap = cmap((64-n+1):64,:);
end