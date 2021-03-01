function h = meeg_diagnostics_source(data,diag,varargin)

STATFIELDS = {'stat' 'mask'};

if nargin >= 3, figtitle = varargin{1}; else, figtitle = 'Sample'; end
if nargin >= 4, savepath = varargin{2}; else, savepath = ''; end

if ~iscell(data), data = {data}; end

plotcfg = keepfields(diag,{'parameter' 'view'});
plotcfg.latency = diag.snapshottwoi./1000; % second

% extract data if needed
for s = 1:numel(data)
    if isfield(data{s},'avg')
        if ~isfield(data{s},plotcfg.parameter), data{s}.(plotcfg.parameter) = data{s}.avg.(plotcfg.parameter); end
        data{s} = rmfield(data{s},'avg');
    end
end

if isfield(data{1},'dim')
    dimord = data{1}.dimord;
    cfg = [];
    cfg.downsample = 2;
    cfg.parameter = plotcfg.parameter;
    for s = 1:numel(data)
        data{s} = ft_sourceinterpolate(cfg, data{s}, diag.mri);
        data{s}.dimord = dimord;
    end
elseif isfield(data{1},'tri') && isfield(diag,'surface')
    dimord = data{1}.dimord;
    cfg = [];
    cfg.parameter = plotcfg.parameter;
    cfg.interpmethod = 'nearest';
    if isfield(data{1},'stat')
        stat = keepfields(data{1}.stat,STATFIELDS);
        data{1} = struct_update(data{1},stat);
    end
    for s = 1:numel(data)
        tmpcfg = cfg;
        if s == 1 && isfield(data{1},'stat'), tmpcfg.parameter = [tmpcfg.parameter STATFIELDS]; end
        data{s} = ft_sourceinterpolate(tmpcfg, data{s}, diag.surface);
        if s == 1 && ~isempty(setdiff(tmpcfg.parameter,cfg.parameter))
            stat = [];
            for f = STATFIELDS, stat.(f{1}) = data{1}.(f{1}); end
            data{1}.stat = stat;
        end
        data{s}.dimord = dimord;
    end
end

fig = meeg_plot(plotcfg,data);

set(fig,'Name',figtitle)

if nargin == 4 && ~isempty(savepath)
    figFn = savepath;
    print(fig,'-noui',[figFn '_multiplot.jpg'],'-djpeg','-r300');
    close(fig);
end
end
