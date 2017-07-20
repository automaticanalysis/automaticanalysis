function [list ll] = dir_wfp(s,dirs)
%%% This is a wrapper around dir and returns filenames plus specified path
%%% prefixes to the files names and returns the output as a cell array.
%%% This is largely just a convenience for quickly getting cell array lists
%%% of file names to pass into the GLM_Flex analysis routines.
%%%
%%%
%%%
%%% Written by Aaron Schultz (aschultz@martinos.org)
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
if nargin==1
    dirs = 0;
end
nn = dir(s);
list = cell(numel(nn),1);
ll = cell(numel(nn),1);
ind = find(s==filesep);
if ~isempty(ind)
    path = s(1:ind(end));
else
    path = [];
end
c = 0;
for ii = 1:length(nn);
    switch dirs
        case 0
            c = c+1;
            list{c,1} = [path nn(ii).name];
            ll{c,1} = nn(ii).name;
        case 1
            if nn(ii).isdir==1
                c = c+1;
                list{c,1} = [path nn(ii).name];
                ll{c,1} = nn(ii).name;
            end
        case -1
            if nn(ii).isdir==0
                c = c+1;
                list{c,1} = [path nn(ii).name];
                ll{c,1} = nn(ii).name;
            end
    end
end
list = list(1:c);
ll = ll(1:c);
