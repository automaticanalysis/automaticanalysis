function dist = DistMat(RS,opt)
%%% Compute the euclidean distance between each point in RS
%%%
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

if nargin == 1
    opt = 1;
end

dist = zeros(length(RS),length(RS));

for ii = 1:size(RS,1);
    TmpD = zeros(size(RS));
    for kk = 1:size(RS,2);
        TmpD(:,kk) = (RS(:,kk)-RS(ii,kk)).^2;
    end
    
    TmpD = sqrt(sum(TmpD,2));
    
    dist(:,ii) = TmpD;
end

if opt == 1
    dist = dist./max(max(dist));
end

if opt == 2
    dist = dist./mean(dist(1:numel(dist)));
end

