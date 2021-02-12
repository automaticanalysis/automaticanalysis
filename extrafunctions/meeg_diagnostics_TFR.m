function h = meeg_diagnostics_TFR(tfr,diag,figtitle,savepath)
TOPO_MAXCOL = 4; % maximum number of columns of topoplot mosaic 
TOL = 1e-8;

if nargin < 3, figtitle = 'Sample'; end

if ~iscell(tfr), tfr = {tfr}; end

if ~isfield(diag,'layout'), diag.layout = ft_prepare_layout([],tfr{1}); end
dispcfg = [];
dispcfg.parameter = diag.parameter;
dispcfg.layout = diag.layout;
dispcfg.marker = 'labels';
dispcfg.interactive = 'no';
dispcfg.showoutline = 'yes';
dispcfg.comment = 'no';
try
    dispcfg.zlim = [0 prctile(cell2mat(cellfun(@(x) x.(diag.parameter), tfr, 'UniformOutput', false)),99,'all')]; % max colour at 99 percentile of all data
catch
    allpow = cell2mat(cellfun(@(x) x.(diag.parameter), tfr, 'UniformOutput', false));
    dispcfg.zlim = [0 prctile(allpow(:),99)]; % max colour at 99 percentile of all data
end
dispcfg.colorbar = 'yes';

cfgmulti = dispcfg;
cfgtopo = dispcfg;

cfgmulti.showlabels  = 'yes';
if isfield(tfr{1},'stat')
    stat = tfr{1}.stat;
    tfr{1} = rmfield(tfr{1},'stat');
    if strcmp(stat.cfg.correctm,'cluster')
        if isfield(stat,'posclusters') && ~isempty(stat.posclusters)
            pos_cluster_pvals = [stat.posclusters(:).prob];
            pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
            pos = ismember(stat.posclusterslabelmat, pos_signif_clust);
        else
            pos = stat.mask*0;
        end
        
        if isfield(stat,'negclusters') && ~isempty(stat.negclusters)
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
    tfr{1}.mask = stat.mask;
    cfgmulti.maskparameter = 'mask';
    cfgmulti.maskstyle = 'opacity';
    cfgtopo.highlight = 'on';
    [iLabER,iLabStat] = match_str(tfr{1}.label, stat.label);
    if isfield(diag,'topohighlight')
        topoHighlightFun = str2func(diag.topohighlight);
    else
        topoHighlightFun = @all;
    end
end

if isfield(tfr{1},'time') && numel(tfr{1}.time) > 1 && numel(tfr{1}.freq) > 1
    ptfr = tfr;
    if numel(ptfr) > 1
        ptfr = {ft_math(struct('parameter',diag.parameter,'operation','subtract'),ptfr{:})};
        if isfield(tfr{1},'mask'), ptfr{1}.mask = tfr{1}.mask; end
        cfgmulti.zlim = prctile(ptfr{1}.(diag.parameter)(:),[1 99]);
    end
    ptfr = ptfr{1};
    if numel(ptfr.label) == 1 % single channel data
        cfgmulti.linewidth = 2;
        ft_singleplotTFR(cfgmulti, ptfr);
        diag.videotwoi =  [];
        diag.snapshottwoi = [];
    else
        ft_multiplotTFR(cfgmulti, ptfr);
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
elseif ~isfield(tfr{1},'time') || numel(tfr{1}.time) == 1
    diag.videotwoi =  [];
    for t = 1:numel(tfr)
        if ~isfield(tfr{t},'time'), tfr{t}.time = 0; end
    end
    diag.snapshottwoi = [tfr{1}.time tfr{1}.time]*1000;
end

if numel(tfr) > 1
    labels = arrayfun(@(x) sprintf('group #%d',x),1:numel(tfr),'UniformOutput',false);
else
    labels = {''};
end

if exist('stat','var')
    tfr = [tfr(1) tfr];
    tfr{1}.(diag.parameter) = stat.stat;
    labels = [{'stat'} labels];
end
    
if ~isempty(diag.snapshottwoi) && ~isempty(diag.snapshotfwoi)
    indTFR = 1:numel(tfr);
    if strcmp(labels{1},'stat'), indTFR = 2:numel(tfr); end
    
    % specify color scaling based on data
    alldat = cellfun(@(x) x.(diag.parameter), tfr(indTFR), 'UniformOutput', false);
    alldat = vertcat(alldat{:});
    if numel(tfr{1}.time) > 1 || numel(tfr{1}.freq) > 1 % smooth multiD data
        if numel(tfr{1}.time) > 1 % along time
            zlimSmooth = mean(diff(diag.snapshottwoi'))/1000*1/mean(diff(tfr{1}.time))*3; % increase the size for the gaussian kernel
        else % along frequencies
            zlimSmooth = mean(diff(diag.snapshotfwoi'))/mean(diff(tfr{1}.freq))*3; % increase the size for the gaussian kernel
        end        
        alldat = arrayfun(@(x) smoothdata(alldat(x,:,:),'gaussian',zlimSmooth),1:size(alldat,1),'UniformOutput',false);
        alldat = horzcat(alldat{:});
    end
    zlim = prctile(alldat(:),[1 99]);
   
    for f = 1:size(diag.snapshotfwoi,1)
        cfgtopo.ylim = diag.snapshotfwoi(f,:);
        fig = figure('Name',sprintf('%s_freq-%1.2f-%1.2f',figtitle,diag.snapshotfwoi(f,:)));
        colPlot = min([TOPO_MAXCOL numel(tfr)*ceil(sqrt(size(diag.snapshottwoi,1)))]);
        colPlot = floor(colPlot/numel(tfr))*numel(tfr);
        rowPlot = ceil(size(diag.snapshottwoi,1)*numel(tfr)/colPlot);
        if rowPlot > colPlot % figure will likely not fit to screen
            warning('Topoplot figure is too large to fit to screen -> It will be hidden.')
            set(fig,'Visible','off');
        end
        set(fig,'Position',[0,0,1080 ceil(1080/colPlot*rowPlot)]);
        
        cbOn = false;
        for t = 1:size(diag.snapshottwoi,1)
            cfgtopo.xlim = diag.snapshottwoi(t,:)/1000;
            
            if exist('stat','var')
                % Get the index of each significant channel
                sampStart = find(tfr{1}.time >= cfgtopo.xlim(1),1,'first');
                sampEnd = find(tfr{1}.time <= cfgtopo.xlim(2),1,'last');
                freqStart = find(stat.freq >= cfgtopo.ylim(1),1,'first');
                freqEnd = find(stat.freq <= cfgtopo.ylim(2),1,'last');
                
                pos_int = zeros(numel(tfr{1}.label),1);
                neg_int = zeros(numel(tfr{1}.label),1);
                pos_int(iLabER) = topoHighlightFun(squeeze(sum(pos(iLabStat, freqStart:freqEnd, sampStart:sampEnd),2)), 2);
                neg_int(iLabER) = topoHighlightFun(squeeze(sum(neg(iLabStat, freqStart:freqEnd, sampStart:sampEnd),2)), 2);
                cfgtopo.highlightchannel = find(pos_int | neg_int);
            end
            
            for e = 1:numel(tfr)
                if all(cfgtopo.xlim < min(tfr{e}.time)) || all(cfgtopo.xlim > max(tfr{e}.time)) ||...
                        all(cfgtopo.ylim < min(tfr{e}.freq)) || all(cfgtopo.ylim > max(tfr{e}.freq))
                    warning('no data within the range')
                    continue
                end
                tmpdat = ft_selectdata(struct('latency',cfgtopo.xlim,'frequency',cfgtopo.ylim),tfr{e});
                if all(isnan(tmpdat.(diag.parameter)(:))), continue; end % skip empty data
                
                cfgtopo.colorbar = 'no';
                if strcmp(labels{e},'stat')
                    cfgtopo.zlim = prctile(tfr{e}.(diag.parameter)(:),[1 99]);
                    cfgtopo.colorbar = 'West';
                else
                    cfgtopo.zlim = zlim;
                    if ~cbOn && (isempty(labels{e}) || strcmp(labels{e},'group #1')), cfgtopo.colorbar = 'West'; cbOn = true; end
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
                ax(t,e) = subplot(rowPlot,colPlot,(t-1)*numel(tfr)+e);
                title(sprintf('%s %03d-%03d ms',labels{e},diag.snapshottwoi(t,:)));
                ft_topoplotTFR(cfgtopo, tfr{e});
            end
        end
        % Set colormaps
        for t = 1:size(diag.snapshottwoi,1)
            for e = 1:numel(tfr)
                if ~isempty(cmaps{t,e}), colormap(ax(t,e),cmaps{t,e}); end
            end
        end        
        % Put colorbars between the axes
        cbs = findall(fig,'type','ColorBar');
        axs = findall(fig,'type','Axes');
        dPos = cbs(1).Position.*[0 0 0.25 0.25];
        if numel(axs) > 1
            dPos(1) = (-(diff(arrayfun(@(x) x.Position(1), axs(1:2))))-axs(1).Position(3))/2;
        end
        for c = cbs'
            c.YAxisLocation = 'left';
            c.Location = 'manual';
            c.Position = c.Position-dPos;
        end
        if nargin >= 4
            set(fig,'PaperPositionMode','auto');
            print(fig,'-noui',sprintf('%s_topoplot_freq-%1.2f-%1.2f.jpg',savepath,diag.snapshotfwoi(f,:)),'-djpeg','-r150');
            close(fig)
        else
            h.topoplot(f) = fig;
        end
    end
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