function groupDesc = ft_granddescriptives(cfg,varargin)

% FT_GRANDDESCRIPTIVES computes the descriptives (mean and SEM) of the specified measure over multiple subjects
%
% Use as
%   [groupDesc] = ft_granddescriptives(cfg, data1, data2, data3...)
%
% The input data data1..N are obtained from either FT_TIMELOCKANALYSIS or FT_FREQANALYSIS and MUST have consistent dimords.
% The configuration structure can contain
%   cfg.channel   = Nx1 cell-array with selection of channels (default = 'all'),
%                   see FT_CHANNELSELECTION for details
%   cfg.latency   = [tmin tmax] or 'all', to specify a subset of latencies (default = 'all')
%   cfg.parameter = string or cell-array of strings indicating which
%                   parameter(s) to average.
%
% See also FT_TIMELOCKGRANDAVERAGE, FT_FREQANALYSIS, FT_FREQDESCRIPTIVES,
% FT_FREQBASELINE

% Copyright (C) 2020, Tibor Auer

%% Configuration
cfg.parameter = ft_getopt(cfg, 'parameter',  []);
if isempty(cfg.parameter)
    ft_error('you should specify a valid parameter to compute');
end
if ischar(cfg.parameter)
    cfg.parameter = {cfg.parameter};
end

groupDesc = keepfields(ft_selectdata(cfg,varargin{1}),...
    {'label','elec','dimord','time','freq','inside','pos','tri',cfg.parameter{1}});
initdata = zeros(size(groupDesc.(cfg.parameter{1})));

%% Computation
for k=1:numel(cfg.parameter)
    outsum = initdata;
    outssq = initdata;
    n      = initdata;
    for j = varargin
        tmp = ft_selectdata(cfg,j{1});
        tmp = tmp.(cfg.parameter{k});
        n      = n + double(isfinite(tmp));
        tmp(~isfinite(tmp)) = 0;
        outsum = outsum + tmp;
        outssq = outssq + tmp.^2;
    end
    
    groupDesc.(cfg.parameter{k}) = outsum./n;
    groupDesc.([cfg.parameter{k} 'sem']) = sqrt((outssq - (outsum.^2)./n)./(n - 1)./n);
end