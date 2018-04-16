function m = crtlFor(m,x,opt)
%%% Written by Aaron Schultz (aschultz@martinos.org)
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
%%% GNU General Public License for more details.

mn = mean(m);
if numel(size(m))==4
    m = reshapeWholeBrain(size(m),m);
end

xx = [ones(size(x,1),1) x];
b = pinv(xx)*m;
% b = (xx\m);

pred = xx*b;
m = m-pred;

if nargin>2
    m = m+repmat(mn,size(m,1),1);
end

% m = m+repmat(tmp,size(m,1),1);