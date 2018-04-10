function [in in2] = matchAll(pl, rv, trm);
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

if nargin == 2;
    trm = [0 0];
end

in = []; 
for ii = 1:length(pl);
    c = 0;
    for jj = 1:length(rv)
        if ~isempty(strfind(rv{jj},pl{ii}(1+trm(1):end-trm(2))))
            c = c+1;
            in(ii,c) = jj;
            in2(jj,c) = ii;
        end
    end
end
