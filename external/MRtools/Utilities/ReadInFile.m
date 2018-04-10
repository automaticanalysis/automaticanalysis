function M = ReadInFile(filename,delim, opt)
%%% Read in a .csv file as a matlab array. If opt = 1 then numbers will be
%%% converted into number format else everything is read in as a string.
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


if nargin == 2;
    opt = 0;
end
fid = fopen(filename);
% try

M1 = textscan(fid, '%s', 'delimiter', '\n'); %, 'bufsize', 16000
% catch
%     M1 = textscan(fid, '%s', 'delimiter', '\n', 'bufsize', 1600000000);
% end
M1 = M1{1};
fclose(fid);

%%%  Read In Data Set from csv
M = [];
for ii = 1:length(M1);
    tmp = M1{ii};
    if delim == ' ';
        ti = regexp(tmp,'  ');
        while numel(ti)>0;
            tmp(ti) = '';
            ti = regexp(tmp,'  ');
        end
    end
    ind = regexp(tmp,delim);
    if isempty(ind)
        M{ii,1} = tmp;
        continue;
    end
    
    ind = [0 ind length(tmp)+1];
    for jj = 1:length(ind)-1;
        
        tmp3 = tmp(ind(jj)+1:ind(jj+1)-1);
        
        if opt == 1;
            if any(tmp3=='/');
                tmp3 = strtrim(tmp3);
            elseif strcmpi(tmp3,'\N');
                tmp3 = NaN;
            elseif isnan(str2num(tmp3)) 
                tmp3 = str2num(tmp3);
            elseif ~isempty(logical(str2num(tmp3)));
                tmp3 = str2num(tmp3);
            else
                tmp3 = strtrim(tmp3);
            end
            if isempty(tmp3)
%                 tmp3 = 'NA';
                tmp3 = NaN;
            end
        end
        M{ii,jj} = tmp3;
    end
end