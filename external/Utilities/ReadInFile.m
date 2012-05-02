function M = ReadInFile(filename,delim, opt)

%%% Read in a .csv file as a matlab array. If opt = 1 then numbers will be
%%% converted into number format else everything is read in as a string.
%%%
%%%  Aaron P. Schultz (ap.schultz@gmail.com) Nov. 10th, 2009


if nargin == 2;
    opt = 0;
end
fid = fopen(filename);
% try
    M1 = textscan(fid, '%s', 'delimiter', '\n', 'bufsize', 16000);
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
        end
        M{ii,jj} = tmp3;
    end
end