function h = meeg_diagnostics_source(source,diag,figtitle,savepath)
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

cfgdiag.funparameter  = 'pow';

if nargin < 3, figtitle = 'Sample'; end

if ~iscell(source), source = {source}; end

for s = 1:numel(source)
    if ~isfield(source{s},'avg')
        source{s}.avg.(cfgdiag.funparameter) = source{s}.(cfgdiag.funparameter);
    end
end

if numel(source) > 1
    labels = arrayfun(@(x) sprintf('group #%d',x),1:numel(source),'UniformOutput',false);
else
    labels = {''};
end

if isfield(source{1},'stat')
    stat = source{1}.stat;
    source{1} = rmfield(source{1},'stat');
%     pos = stat.mask & (stat.stat > 0);
%     neg = stat.mask & (stat.stat < 0);
    source = [source(1) source];
    source{1}.avg.(cfgdiag.funparameter) = stat.stat;
    for s = 1:numel(source)
        source{s}.mask = stat.mask;
    end
    TOPO_MAXCOL = TOPO_MAXCOL/2;
    labels = [{'stat'} labels];
end

if isfield(source{1},'dim'), TOPO_MAXCOL = TOPO_MAXCOL/2; end
if isfield(source{1},'time')
    if isempty(diag.snapshottwoi), diag.snapshottwoi = [source{1}.time(1) source{1}.time(end)]; end
    sizeplot = [ceil(size(diag.snapshottwoi,1)/TOPO_MAXCOL) min(size(diag.snapshottwoi,1),TOPO_MAXCOL)*numel(source)];
else
    diag.snapshottwoi = [0 0];
    for s = 1:numel(source)
        source{s}.time = 0;
    end
    sizeplot = [1,numel(source)];
end

if isfield(source{1},'dim')
    sizeplot(2) = sizeplot(2)*2;
    cfg = [];
    cfg.downsample = 2;
    cfg.parameter = cfgdiag.funparameter;
    for s = 1:numel(source)
        source{s} = ft_sourceinterpolate(cfg, source{s}, diag.mri);
    end
end

for s = 1:numel(source)
    minval(s) = min(source{s}.avg.(cfgdiag.funparameter)(:));
    maxval(s) = max(source{s}.avg.(cfgdiag.funparameter)(:));
    % colormap for surface
    if (minval(s) < 0) && (maxval(s) > 0), cmap{s} = [winter; hot];
    elseif minval(s) < 0, cmap{s} = winter;
    else, cmap{s} = hot;
    end
end

for f = 1:size(diag.snapshotfwoi,1)
    fig = figure;
    set(fig, 'color', [1 1 1]);
    for t = 1:size(diag.snapshottwoi,1)
        for s = 1:numel(source)
            p = subplot(sizeplot(1),sizeplot(2),(t-1)*numel(source)+s);
            strTitle = sprintf('%03d-%03d ms',diag.snapshottwoi(t,:));
            if ~isempty(labels{s}), strTitle = sprintf('%s: %s',labels{s},strTitle); end
            title(strTitle);
            cfg = [];
            cfg.frequency = diag.snapshotfwoi(f,:);
            cfg.avgoverfreq = 'yes';
            cfg.latency = diag.snapshottwoi(t,:)/1000;
            cfg.avgovertime = 'yes';
            sourcePlot = ft_selectdata(cfg,source{s});
            sourcePlot.avg.(cfgdiag.funparameter) = sourcePlot.(cfgdiag.funparameter);
            
            if isfield(sourcePlot,'dim')
                %             p = subplot(sizeplot(1),sizeplot(2),(t-1)*2+1);
                %             cfg = cfgdiag;
                %             cfg.title = 'positive';
                %             cfg.maskparameter = cfg.funparameter;
                %             cfg.funcolorlim   = [0 maxval];
                %             cfg.opacitylim    = cfg.funcolorlim ;
                %             cfg.method = 'slice';
                %             ft_sourceplot(cfg, sourcePlot, diag.mri);
                %             copyobj(gca,p);
                %
                %             p = subplot(sizeplot(1),sizeplot(2),(t-1)*2+2);
                %             cfg = cfgdiag;
                %             cfg.title = 'negative';
                %             cfg.maskparameter = cfg.funparameter;
                %             cfg.funcolorlim   = [minval 0];
                %             cfg.opacitylim    = cfg.funcolorlim ;
                %             cfg.method = 'slice';
                %             ft_sourceplot(cfg, sourcePlot, diag.mri);
                %             arrayfun(@(x) copyobj(x,p), get(gca,'Children'));
                
            elseif isfield(sourcePlot,'tri')
                if ~isfield(sourcePlot,'mask')
                    mask = abs(sourcePlot.avg.(cfgdiag.funparameter)); mask = mask./max(mask);
                    sourcePlot.mask = mask;
                else % stats
                    if ~any(sourcePlot.mask(:)), title(sprintf('%s\nns.',strTitle)); end
                end
                
                cfg = cfgdiag;
                cfg.funcolormap = cmap{s};
                cfg.funcolorlim = [minval(s) maxval(s)];
                cfg.maskparameter = 'mask';
                cfg.opacitylim = cfg.funcolorlim;
                cfg.method = 'surface';
                cfg.camlight = 'no';
                ft_sourceplot(cfg, sourcePlot);
                arrayfun(@(x) copyobj(x,p), get(gca,'Children'));
                set(p,'CLim',get(gca,'CLim'))
                view(p,VAs{strcmp(VAs(:,1),diag.view),2}{:});
                camlight(p);
                axis(p,'image');
                set(p, 'Tag', 'ik', 'Visible', 0);
                set(p, 'Tag', 'jk', 'Visible', 0);
                set(p, 'Tag', 'ij', 'Visible', 0);
                set(get(p,'Title'), 'Visible', 1);
                colormap(p,colormap(gca));
                colorbar(p);
                close(gcf);
            end
        end
    end
    set(fig,'Name',sprintf('%s_freq-%1.2f-%1.2f',figtitle,diag.snapshotfwoi(f,:)));
    
    if nargin >= 4
        set(fig,'Position',[0,0,1080 1080]);
        set(fig,'PaperPositionMode','auto');
        print(fig,'-noui',sprintf('%s_source_freq-%1.2f-%1.2f.jpg',savepath,diag.snapshotfwoi(f,:)),'-djpeg','-r300');
        close(fig);
    else
        h(f) = fig;
    end
end
