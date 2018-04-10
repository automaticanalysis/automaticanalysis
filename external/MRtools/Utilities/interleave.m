function out = interleave(varargin)
%%% A function to interleave a set of numeric values or cells.
%%% This function gets used to interleave filenames to get the correct
%%% order of specification for GLM_Flex wrappers
%%%
%%% Inputs
%%% a separate entry for each group to interleave
%%%
%%% Examples:
%%%
%%% interleave([1 3 5], [2 4 6])
%%%     produces [1 2 3 4 5 6]'
%%% 
%%% interleave({'a' 'b' 'c'}, {'d' 'e' 'f'}, {'g' 'h' 'i'})
%%%     produces {'a' 'd' 'g' 'b' 'e' 'h' 'c' 'f' 'i'}
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


input = [];
if numel(varargin)==1
    input = varargin{1};
else
    for ii = 1:numel(varargin)
        input{ii} = varargin{ii};
    end
end

lens = [];
for ii = 1:length(input)
    lens(ii) = numel(input{ii});
end

if mean(diff(lens))~=0;
    out = [];
    return
end

if iscell(input{1})
    out = cell(sum(lens),1);
else
    out = zeros(sum(lens),1);
end

for ii = 1:length(input);
    out(ii:length(input):numel(out)) = input{ii};
end