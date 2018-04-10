function [cols cm cc] = cmap(X, lims, cm, nBins)
%%% Written by Aaron P. Schultz - aschultz@martinos.org
%%%
%%% Copyright (C) 2014,  Aaron P. Schultz
%%%
%%% Supported in part by the NIH funded Harvard Aging Brain Study (P01AG036694) and NIH R01-AG027435 
%%%
%%% This program is free software: you can redistribute it and/or modify
%%% it under the terms of the GNU General Public License as published by
%%% the Free Software Foundation, either version 3 of the License, or
%%% any later version.
%%% 
%%% This program is distributed in the hope that it will be useful,
%%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%% GNU General Public License for more details.X = X(:);

lims = sort(lims);
if nargin<4
    nBins = 256;
end

if ischar(cm)
    eval(['cm = colmap(''' cm ''',' num2str(nBins) ');']);
end

nBins = size(cm,1);

X((X<lims(1)))=lims(1);
X((X>lims(2)))=lims(2);


cc = [X(:); lims(:)];
cc = cc/(diff(lims));
cc = (cc*(nBins));
cc = cc+(nBins-max(cc));
cc = floor(cc(1:end-2))+1;
cc(cc>nBins)=nBins;
cc(cc<1)=1;

cols = nan(numel(cc),3);
cols(~isnan(cc),:) = cm(cc(~isnan(cc)),:);

