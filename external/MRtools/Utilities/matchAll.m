function [in ord] = matchAll(pl, rv, trm);
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
if nargin == 2;
    trm = [0 0];
end

in = [];
ord = [];

c = 0;
for ii = 1:length(pl);
    %tmp = strmatch(pl{ii}(1+trm(1):end-trm(2)),rv);
    tmp = searchCellStr(['^' pl{ii}(1+trm(1):end-trm(2))  '$'],rv);
    if ~isempty(tmp)
        c = c+1;
        in(c) = ii; 
        try
            ord(c) = tmp;
        catch
            warning(['duplicate found for ' rv{tmp(1)} ' using first instance']);
            ord(c) = tmp(1);
        end
    end
end