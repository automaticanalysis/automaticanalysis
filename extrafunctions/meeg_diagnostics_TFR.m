function fig = meeg_diagnostics_TFR(data,diag,varargin)

if nargin >= 3, figtitle = varargin{1}; else, figtitle = 'Sample'; end
if nargin >= 4, savepath = varargin{2}; else, savepath = ''; end

if ~iscell(data), data = {data}; end

plotcfg.parameter = 'powspctrm';
plotcfg.latency = diag.snapshottwoi./1000; % second

% tweak the data into TFR
if numel(data{1}.label) < 4 % plot single channels
    plotcfg.channels = unique(data{1}.label);
end
plotcfg.layout = ft_prepare_layout([],data{1});

% plot as TFR
fig = meeg_plot(plotcfg,data);

set(fig,'Name',figtitle)

if nargin == 4 && ~isempty(savepath)
    figFn = savepath;
    if ~contains(figFn,'plot'), figFn = [figFn '_multiplot.jpg']; end
    print(fig,'-noui',figFn,'-djpeg','-r300');
    close(fig);
end

if isfield(diag,'snapshotfwoi') && ~isempty(diag.snapshotfwoi) && isfield(data{1},'freq') && numel(data{1}.freq) > 1
    for f = 1:size(diag.snapshotfwoi,1)
        for d = 1:numel(data)
            dat{d} = ft_selectdata(struct('frequency',diag.snapshotfwoi(f,:),'avgoverfreq','yes'),data{d});
            dat{d}.dimord = strrep(dat{d}.dimord,'_freq','');
        end
        fig = meeg_plot(plotcfg,dat);
        
        set(fig,'Name',figtitle)

        if nargin == 4 && ~isempty(savepath)
            figFn = savepath;
            print(fig,'-noui',[figFn '_topoplot_freq-' sprintf('%1.2f-%1.2f',diag.snapshotfwoi(f,:)) '.jpg'],'-djpeg','-r300');
            close(fig);
        end
    end
end

end

