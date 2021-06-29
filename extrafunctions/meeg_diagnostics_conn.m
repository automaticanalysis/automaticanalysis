function fig = meeg_diagnostics_conn(data,diag,varargin)

if nargin >= 3, figtitle = varargin{1}; else, figtitle = 'Sample'; end
if nargin >= 4, savepath = varargin{2}; else, savepath = ''; end
if ~isfield(diag,'parameter'), diag.parameter = 'crsspctrm'; end

if ~iscell(data), data = {data}; end

if (ndims(data{1}.(diag.parameter)) > 3 && ~any(strcmp(strsplit(data{1}.dimord,'_'),'time'))) || ndims(data{1}.(diag.parameter)) > 4
    aas_log([],false,'WARNING: visualising data with more than 3D+time is not yet implemented')
    return
end

plotcfg.parameter = diag.parameter;
plotcfg.latency = diag.snapshottwoi./1000; % second

% tweak the data into freq data
cf = cellfun(@cf2tf, data,'UniformOutput',false);

% plot as TFR
for p = 1:numel(cf{1})
    tf = cellfun(@(x) x(p), cf, 'UniformOutput',false);

    fig(p) = meeg_plot(plotcfg,tf);    
    
    set(fig(p),'Name',figtitle)
    
    if nargin == 4 && ~isempty(savepath)
        figFn = savepath;
        print(fig(p),'-noui',[figFn '_multiplot.jpg'],'-djpeg','-r300');
        close(fig(p));
    end
end
end

function tf = cf2tf(data)
TFDIMS = {'freq' 'rpt' 'rpttap'};
dims = strsplit(data.dimord,'_'); 
dims(1) = []; % not for chancmb
dims(strcmp(dims,'time')) = []; % not for time, if exists
tf = data;
tf.dimord = 'chancmb';
for d = 1:numel(dims)
    tf.(TFDIMS{d}) = data.(dims{d});
    tf.dimord = strjoin([tf.dimord TFDIMS(d)],'_');
end
tf = rmfield(tf,setdiff(dims,TFDIMS));
if isfield(tf,'time'), tf.dimord = strjoin({tf.dimord 'time'},'_'); end
end
