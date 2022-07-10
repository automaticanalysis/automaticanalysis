function h = meeg_diagnostics_source(data,diag,varargin)

STATFIELDS = {'stat' 'mask'};

if nargin >= 3, figtitle = varargin{1}; else, figtitle = 'Sample'; end
if nargin >= 4, savepath = varargin{2}; else, savepath = ''; end

if ~iscell(data), data = {data}; end

plotcfg = keepfields(diag,intersect(fieldnames(diag),{'parameter' 'view' 'statlim'}));
plotcfg.latency = diag.snapshottwoi./1000; % second

% extract data if needed
for s = 1:numel(data)
    if isfield(data{s},'avg')
        if ~isfield(data{s},plotcfg.parameter), data{s}.(plotcfg.parameter) = data{s}.avg.(plotcfg.parameter); end
        data{s} = rmfield(data{s},{'avg' 'method'});
    end
    data{s} = rmfield(data{s},intersect(fieldnames(data{s}),{'aparc','aparclabel'}));
end

if isfield(diag,'background')
    dimord = data{1}.dimord;
    cfg = keepfields(plotcfg,'parameter');
    if isfield(data{1},'stat')
        stat = keepfields(data{1}.stat,STATFIELDS);
        data{1} = struct_update(data{1},stat);
    end
    
    if isfield(data{1},'dim')
        cfg.interpmethod = 'linear';
        cfg.downsample = 2;
    elseif isfield(data{1},'tri')
        cfg.interpmethod = 'nearest';
    end
    
    for s = 1:numel(data)
        tmpcfg = cfg;
        if s == 1 && isfield(data{1},'stat'), tmpcfg.parameter = [tmpcfg.parameter STATFIELDS]; end
        data{s} = ft_sourceinterpolate(tmpcfg, data{s}, diag.background);
        if s == 1 && ~isempty(setdiff(tmpcfg.parameter,cfg.parameter))
            stat = [];
            for f = STATFIELDS, stat.(f{1}) = data{1}.(f{1}); end
            stat.mask(isnan(stat.mask)) = 0;
            data{1}.stat = stat;
            data{1} = rmfield(data{1},'mask');
        end
        data{s}.dimord = dimord;
    end
end

fig = meeg_plot(plotcfg,data);

set(fig,'Name',figtitle)

if nargin == 4 && ~isempty(savepath)
    figFn = savepath;
    print(fig,'-noui',[figFn '.jpg'],'-djpeg','-r300');
    close(fig);
end
end
