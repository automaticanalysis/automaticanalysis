function out = ft_combine(cfg,varargin)
% FT_COMBINE combines FieldTrip data structures based on weights
%
% You can specify the following configuration options:
%   cfg.parameter = string indicating which parameter to analyse. default is set to 'avg', if it is present in the data
%   cfg.weights   = 1xN array of weights, where N = number of inputs (default = ones(1,N))
%   cfg.normalise = 'no' or 'yes'  normalise weights (default = 'no')

% Copyright (C) 2020, Tibor Auer
%
% $Id$

% set the defaults
cfg.parameter = ft_getopt(cfg, 'parameter', 'avg');
cfg.weights = ft_getopt(cfg, 'weights', ones(1,numel(varargin)));
cfg.normalise = ft_getopt(cfg, 'normalise', 'no');

if strcmp(cfg.normalise,'yes')
    cfg.weights(cfg.weights>0) = cfg.weights(cfg.weights>0)/sum(cfg.weights(cfg.weights>0));
    cfg.weights(cfg.weights<0) = -cfg.weights(cfg.weights<0)/sum(cfg.weights(cfg.weights<0));
end
weights = cfg.weights;

data = varargin;
cfg = keepfields(cfg,'parameter'); cfg.operation = 'multiply';
for i = 1:numel(data) 
    cfg.scalar = weights(i);
    out{i} = ft_math(cfg,data{i});
end

if numel(out) > 1
    cfg = keepfields(cfg,'parameter'); cfg.operation = 'add';
    out = ft_math(cfg,out{:});
else
    out = out{1};
end
end