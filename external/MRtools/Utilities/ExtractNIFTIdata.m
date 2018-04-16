function X = ExtractNIFTIdata(fn,ml)
%%% Written by Aaron Schultz (aschultz@martinos.org)
%%%
%%% Last Updated Dec. 11 2012;
%%%
%%% Copyright (C) 2012,  Aaron P. Schultz
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
%%% GNU General Public License for more details.

if ischar(fn)
    h = spm_vol(fn);
end

if iscell(fn)
    h = spm_vol(char(fn));
end

if isstruct(fn)
    h = fn;
end

% [ml vi] = getMatCoord(h(1),[0 -53 26],12);

X = zeros(numel(h),size(ml,1));
for ii = 1:numel(h)
    X(ii,:) = spm_sample_vol(h(ii),ml(:,1),ml(:,2),ml(:,3),0);
end

