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

% band -> freq
if isfield(data{1},'band') % assume same data type
    if isnumeric(diag.snapshotfwoi{1}) % convert band to freq
        for d = 1:numel(data)
            if isfield(data{d},'bandspec')
                bands = data{d}.bandspec.bandbound;
                freq = cellfun(@mean, bands)';
            else
                ft_error('no band specification found, cannot convert band to Hz');
            end
            data{d} = rmfield(data{d},intersect(fieldnames(data{d}),{'band','bandspec'}));
            data{d}.freq = freq;
            data{d}.dimord = strrep(data{d}.dimord,'band','freq');
        end
    elseif ischar(diag.snapshotfwoi{1})
        for d = 1:numel(data)
            [~,indDiag,indData] = intersect(diag.snapshotfwoi,data{1}.band,'stable');
            data{d} = rmfield(data{d},intersect(fieldnames(data{d}),{'band','bandspec'}));
            data{d}.freq = indData;
            data{d}.dimord = strrep(data{d}.dimord,'band','freq');
        end
        diagBands = diag.snapshotfwoi;
        diag.snapshotfwoi = [indDiag-0.25 indDiag+0.25];
    end
end

% plot as TFR
fig = meeg_plot(plotcfg,data);

set(fig,'Name',figtitle)

if nargin == 4 && ~isempty(savepath)
    figFn = savepath;
    if contains(figFn,'plot'), figFn = [figFn '.jpg'];
    else, figFn = [figFn '_multiplot.jpg']; end
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
            if exist('diagBands','var') % bands
                freqspec = diagBands{round(mean(diag.snapshotfwoi(f,:)))};
            else
                freqspec = sprintf('%1.2f-%1.2f',diag.snapshotfwoi(f,:));
            end
            print(fig,'-noui',[figFn '_topoplot_freq-' freqspec '.jpg'],'-djpeg','-r300');
            close(fig);
        end
    end
end

end

