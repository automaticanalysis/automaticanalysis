function M = reshapeWholeBrain(ss,D)
%%% Shape and reshape 4-D matrix to 2-D matrix and then back again.
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

if numel(size(D))==2
    M = reshape(D', ss(1), ss(2), ss(3), ss(4)); 
end

if numel(size(D))==4
     M = double(reshape(D, prod(ss(1:3)),ss(4))');
end
            
            
           