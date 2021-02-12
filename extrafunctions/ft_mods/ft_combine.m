function out = ft_combine(cfg,varargin)
% FT_COMBINE combines FieldTrip data structures based on weights
%
% You can specify the following configuration options:
%   cfg.parameter = string or cell-array indicating which parameter(s) to compute. default is set to 'avg', if it is present in the data
%   cfg.weights   = 1xN array of weights, where N = number of inputs (default = ones(1,N))
%   cfg.normalise = 'no' or 'yes'  normalise weights (default = 'no')

% Copyright (C) 2020, Tibor Auer
%
% $Id$

% set the defaults
cfg.parameter = ft_getopt(cfg, 'parameter', 'avg');
cfg.weights = ft_getopt(cfg, 'weights', ones(1,numel(varargin)));
cfg.normalise = ft_getopt(cfg, 'normalise', 'no');

if ~iscell(cfg.parameter)
  cfg.parameter = {cfg.parameter};
end
if strcmp(cfg.normalise,'yes')
    cfg.weights(cfg.weights>0) = cfg.weights(cfg.weights>0)/sum(cfg.weights(cfg.weights>0));
    cfg.weights(cfg.weights<0) = -cfg.weights(cfg.weights<0)/sum(cfg.weights(cfg.weights<0));
end
weights = cfg.weights;

data = varargin;
out = data{1};
for p = cfg.parameter
    dat = arrayfun(@(x) data{x}.(p{1})*weights(x),1:numel(data),'UniformOutput',false);
    out.(p{1}) = sum(cat(ndims(dat{1})+1,dat{:}),ndims(dat{1})+1);
end
end