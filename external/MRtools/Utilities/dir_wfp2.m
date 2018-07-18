function [list ll] = dir_wfp2(s)
%%% This is a wrapper around dir and returns filenames plus the full path
%%% to the files names and returns the output as a cell array.
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
%s = 'Rest*/ffBPS_motRes_*.nii'

[a b] = system(['ls -1 ' s]);

list = regexp(b,'\n','split');
list = list(1:end-1)';

if strcmpi('ls: No match.',list{1})
    list = [];
end

if nargout == 2
    ll = cell(size(list));
    for ii = 1:numel(list);
        ind = find(list{ii}==filesep);
        ll{ii} = list{ii}(ind(end)+1:end);
    end
end