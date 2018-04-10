function [set other] = InAll(probe,set,other)
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

for ii = 1:numel(set);
    [a b] = matchAll(probe,set{ii});
    probe = probe(a);
end


if nargin < 3
    other = {};
end

for ii = 1:numel(set);
    [a b] = matchAll(probe,set{ii});
    set{ii} = set{ii}(b);
    if nargin == 3
       other{ii} = other{ii}(b);
    else
        other{ii} = b;
    end
end
   

% set
% [char(set{1}) repmat('    ',numel(set{1}),1) char(set{3})]