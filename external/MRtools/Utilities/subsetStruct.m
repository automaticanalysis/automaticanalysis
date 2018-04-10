function S = subsetStruct(S,index)
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


% if ~isstruct(S)
%     warning off
%     tmp = struct(S);
%     SS = [];
%     for ii = 1:numel(tmp.varnames);
%         SS.(tmp.varnames{ii}) = tmp.data{ii};
%     end
%     S = SS;
%     warning on
% end

list = fields(S);

for ii = 1:numel(list);
    S.(list{ii}) = S.(list{ii})(index);
end