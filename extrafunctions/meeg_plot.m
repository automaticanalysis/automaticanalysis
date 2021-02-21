function h = meeg_plot(cfg,data)
% cfg
%   - parameter  - char
%   - latency       - [Nx2], second

%% Contants
FIGWIDTH = 1080;
TOPO_MAXCOL = 4; % maximum number of columns of plot mosaic, should be even

VAs = {...
    'RAS' {135,45};...
    'LAS' {-135,45};...
    'RPS' {45,45};...
    'LPS' {-45,45};...
    'RAI' {135,-45};...
    'LAI' {-135,-45};...
    'RPI' {45,-45};...
    'LPI' {-45,-45};...
    };

%% Groups
if numel(data) > 1
    labels = arrayfun(@(x) sprintf('group #%d',x),1:numel(data),'UniformOutput',false);
else
    labels = {'group'};
end

%% Stats
if isfield(data{1},'stat')
    stat = data{1}.stat;
    data{1} = rmfield(data{1},'stat');
    pos = stat.mask & (stat.stat > 0);
    neg = stat.mask & (stat.stat < 0);
    data = [data(1) data];
    data{1}.(cfg.parameter) = stat.stat;
    for s = 1:numel(data)
        data{s}.mask = logical(stat.mask);
    end
    TOPO_MAXCOL = TOPO_MAXCOL - 1;
    labels = [{'stat'} labels];
end

%% Time
if isfield(data{1},'time')
    if isempty(cfg.latency), cfg.latency = [data{1}.time(1) data{1}.time(end)]; end
    sizeplot = [ceil(size(cfg.latency,1)/TOPO_MAXCOL) min(size(cfg.latency,1),TOPO_MAXCOL)*numel(data)];
else
    cfg.latency = [0 0];
    for s = 1:numel(data)
        data{s}.time = 0;
    end
    sizeplot = [1,numel(data)];
end

%% Scale
for s = 1:numel(data)
    if isfield(data{s},'mask') && any(data{s}.mask(:))
        mask = data{s}.mask;
    else
        mask = true(size(data{s}.(cfg.parameter)));
    end
    if strcmp(labels{s},'stat')
        minval(s) = prctile(data{s}.(cfg.parameter)(mask),1,'all');
        maxval(s) = prctile(data{s}.(cfg.parameter)(mask),99,'all');
    else
        minval(s) = prctile(cell2mat(cellfun(@(d) d.(cfg.parameter)(mask), data(~strcmp(labels,'stat')), 'UniformOutput',false)),1,'all');
        maxval(s) = prctile(cell2mat(cellfun(@(d) d.(cfg.parameter)(mask), data(~strcmp(labels,'stat')), 'UniformOutput',false)),99,'all');
    end
    % colormaps
    if (minval(s) < 0) && (maxval(s) > 0)
        r = maxval(s)/-minval(s);
        if r > 1, cmaps{s} = [cmapcold(round(64/r)); cmaphot(64)];
        else, cmaps{s} = [cmapcold(64); cmaphot(round(r*64))];
        end
    elseif minval(s) < 0, cmaps{s} = cmapcold(64);
    else
        cmaps{s} = cmaphot(64);
    end
end

%% Plot
h = figure;
set(h, 'color', [1 1 1]);
for t = 1:size(cfg.latency,1)
    for s = 1:numel(data)
        p = subplot(sizeplot(1),sizeplot(2),(t-1)*numel(data)+s, 'Parent',h);
        
        tmpcfg = [];
        tmpcfg.latency = cfg.latency(t,:);
        if isnumeric(cfg.latency)
            strTitle = sprintf('%03d-%03d ms',cfg.latency(t,:)*1000);
            tmpcfg.avgovertime = 'yes';
        else % all
            strTitle = 'all';
        end
        if ~isempty(labels{s}), strTitle = sprintf('%s: %s',labels{s},strTitle); end
        dataPlot = ft_selectdata(tmpcfg,data{s});
        
        if isfield(dataPlot,'mask') % stats
            if ~any(dataPlot.mask(:))
                strTitle = sprintf('%s\nns.',strTitle);
            else
                dataPlot.(cfg.parameter) = dataPlot.(cfg.parameter).*dataPlot.mask;
            end
        end
        
        tmpcfg = [];
        tmpcfg.parameter = cfg.parameter;
        if isfield(dataPlot,'elec')
            tmpcfg.zlim = [minval(s) maxval(s)];
            if ft_datatype(dataPlot,'freq') % TFR
                if any(strcmp(strsplit(dataPlot.dimord,'_'),'time')) % multiplot
                    FIGWIDTH = 2*1080;
                    
                    tmpcfg.showlabels = 'yes';
                    tmpcfg.showoutline = 'yes';
                    tmpcfg.interactive = 'no';
                    ft_multiplotTFR(tmpcfg,dataPlot);
                    arrayfun(@(x) copyobj(x,p), get(gca,'Children'));                    
                else %topoplot
                end
            elseif ft_datatype(data,'timelock')  % ER
            end
            
        elseif isfield(dataPlot,'dim')
            %             p = subplot(sizeplot(1),sizeplot(2),(t-1)*2+1);
            %             cfg = cfgdiag;
            %             cfg.title = 'positive';
            %             cfg.maskparameter = cfg.parameter;
            %             cfg.funcolorlim   = [0 maxval];
            %             cfg.opacitylim    = cfg.funcolorlim ;
            %             cfg.method = 'slice';
            %             ft_dataplot(cfg, dataPlot, cfg.mri);
            %             copyobj(gca,p);
            %
            %             p = subplot(sizeplot(1),sizeplot(2),(t-1)*2+2);
            %             cfg = cfgdiag;
            %             cfg.title = 'negative';
            %             cfg.maskparameter = cfg.parameter;
            %             cfg.funcolorlim   = [minval 0];
            %             cfg.opacitylim    = cfg.funcolorlim ;
            %             cfg.method = 'slice';
            %             ft_dataplot(cfg, dataPlot, cfg.mri);
            %             arrayfun(@(x) copyobj(x,p), get(gca,'Children'));
            
        elseif isfield(dataPlot,'tri')
            tmpcfg = cfgdiag;
            tmpcfg.funcolormap = cmaps{s};
            tmpcfg.funcolorlim = [minval(s) maxval(s)];
            tmpcfg.method = 'surface';
            tmpcfg.camlight = 'no';
            ft_sourceplot(tmpcfg, dataPlot);
            arrayfun(@(x) copyobj(x,p), get(gca,'Children'));
            view(p,VAs{strcmp(VAs(:,1),cfg.view),2}{:});
            camlight(p);
        end
        set(p,'CLim',get(gca,'CLim'));
        axis(p,'image');
        set(p, 'Tag', 'ik', 'Visible', 0);
        set(p, 'Tag', 'jk', 'Visible', 0);
        set(p, 'Tag', 'ij', 'Visible', 0);
        title(p, strTitle);
        set(get(p,'Title'), 'Visible', 1);
        colormap(p,cmaps{s});
        colorbar(p);
        close(gcf);
    end
end
set(h,'Position',[0,0,FIGWIDTH,round(sizeplot(1)/sizeplot(2)*FIGWIDTH)]);
set(h,'PaperPositionMode','auto');

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