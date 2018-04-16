function neigh = neighbors(currloc,ss)
%%% Written by Aaron Schultz (aschultz@martinos.org)
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

% currloc = [44 30 16];
% ss = [ 53 63 54];

neigh = zeros(27,3);
cc = 0;
for ii = currloc(1)-1:currloc(1)+1
    for jj = currloc(2)-1:currloc(2)+1
        for kk = currloc(3)-1:currloc(3)+1
            cc = cc+1;
            neigh(cc,1:3) = [ii jj kk];
        end
    end
end