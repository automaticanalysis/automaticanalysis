function dat = GrubbsMat(dat,alpha)
%%% Find outliers using Grubb's test.  This function opperated on columns.
%%%
%%% Input = a data matrix.
%%%
%%% Output = a data matrix where detected outliers are replaced with NaN
%%% values
%%%
%%% Written by Aaron Schultz, May 5th, 2010 (aschultz@martinos.org)
%%% Copyright (C) 2011,  Aaron P. Schultz
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

c = 0;
outind = [];
breaker = 1;
if nargin == 1
    alpha = .05;
end
tail = 2;
c = 0;
while breaker == 1;
    c = c+1;
    mu = nanmean(dat);
    s = nanstd(dat);

    mm = abs(dat-(repmat(mu,size(dat,1),1)));
    [Ma Ii] = max(mm);
    ind = (0:size(dat,1):prod(size(dat))-1) + Ii;

    G = Ma./s;
    N = size(dat,1);

    tcrit=abs(tinv(alpha/(2*N),N-2)); %critical value
    crit = ((N-1)/sqrt(N)) .* sqrt((tcrit.^2)./(N-2+(tcrit.^2)));

    ind2 = find(G>crit);
    if isempty(ind2)
        breaker = 0;
    else
        dat(ind(ind2)) = NaN;
    end
end
