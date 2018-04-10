function [out user] = UserTime
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

tmp = pwd;
cd ~
user = pwd;
cd(tmp);

ind = find(user == filesep);
if ind(end)==numel(user);
    user = user(ind(end-1)+1:ind(end)-1);
else
    user = user(ind(end)+1:end);
end
out = ['Last run by ' user ' on ' datestr(clock)];