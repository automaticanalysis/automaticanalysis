function out = ft_combine(cfg,varargin)
% FT_COMBINE combines FieldTrip data structures based on weights
%
% You can specify the following configuration options:
%   cfg.parameter = string or cell-array indicating which parameter(s) to compute. default is set to 'avg', if it is present in the data
%   cfg.weights   = 1xN array of weights, where N = number of inputs (default = ones(1,N))
%   cfg.contrast  = how to calculate contrast (default = 'difference')
%                       'difference' - works with any number of inputs
%                       'normaliseddifference' - works with any number (>1) of inputs with at least one positive and negative weights
%                       'ratio' - works with only two inputs with negative weight on only one of them
%   cfg.normalise = 'no' or 'yes', normalise weights (default = 'no')

% Copyright (C) 2020, Tibor Auer
%
% $Id$

% set the defaults
cfg.parameter = ft_getopt(cfg, 'parameter', 'avg');
cfg.weights = ft_getopt(cfg, 'weights', ones(1,numel(varargin)));
cfg.contrast = ft_getopt(cfg, 'contrast', 'difference');
if strcmp(cfg.contrast,'normaliseddifference') &&  numel(cfg.weights) > 2, ft_error('normaliseddifference is not implemented for more than two items'); end

cfg.normalise = ft_getopt(cfg, 'normalise', 'no');

if ~iscell(cfg.parameter)
    cfg.parameter = {cfg.parameter};
end
if strcmp(cfg.normalise,'yes')
    cfg.weights(cfg.weights>0) = cfg.weights(cfg.weights>0)/sum(cfg.weights(cfg.weights>0));
    cfg.weights(cfg.weights<0) = -cfg.weights(cfg.weights<0)/sum(cfg.weights(cfg.weights<0));
end
weights = cfg.weights;

% check data
data = varargin;
if numel(data) > 1
    dims = strsplit(data{1}.dimord,'_');
    discrepantDim = logical(diff(cell2mat(cellfun(@(d) size(d.(cfg.parameter{1}))', data, 'UniformOutput', false))'));
    for dim = dims(discrepantDim)
        indDim = find(strcmp(dims,dim{1}));
        scale = cellfun(@(d) d.(dim{1}), data, 'UniformOutput', false);
        sc = scale{1}; for i = 2:numel(scale), sc = intersect(sc, scale{i}); end; scale = sc;        
        for i = 1:numel(data)
            [~,ind] = intersect(data{i}.(dim{1}), scale);
            for p = cfg.parameter
                [dat, shperm, shn] = shiftdata(data{i}.(p{1}),indDim);
                data{i}.(p{1}) = unshiftdata(dat(ind,:,:,:,:),shperm,shn);
            end
            data{i}.(dim{1}) = scale;
        end
    end
end

out = data{1};
for p = cfg.parameter
    switch cfg.contrast
        case {'difference' 'normaliseddifference'}
            dat = arrayfun(@(x) data{x}.(p{1})*weights(x),1:numel(data),'UniformOutput',false);
            out.(p{1}) = sum(cat(ndims(dat{1})+1,dat{:}),ndims(dat{1})+1);
            if strcmp(cfg.contrast,'normaliseddifference')
                out.(p{1}) = out.(p{1})./(-sum(cat(ndims(dat{1})+1,dat{weights<0}),ndims(dat{1})+1));
            end
        case 'ratio'
            dat = arrayfun(@(x) data{x}.(p{1})*abs(weights(x)),1:numel(data),'UniformOutput',false);
            datexp = arrayfun(@(x) ones(size(data{x}.(p{1})))*(weights(x)./abs(weights(x))),1:numel(data),'UniformOutput',false);
            out.(p{1}) = prod(cat(ndims(dat{1})+1,dat{:}).^cat(ndims(datexp{1})+1,datexp{:}),ndims(dat{1})+1);
    end
end
end