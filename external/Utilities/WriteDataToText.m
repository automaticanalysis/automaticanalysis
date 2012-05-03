function WriteDataToText(S, fn, type, sep, pad, percision)
%%%  Handy function  for writing data out to ascii files.
%%%
%%% S is a data structure where each field is a column and each level of
%%% the structure is a row (e.g S(1).col1 = 1; S(2).col1 = 2; S(2).col2 = 2 etc.)
%%%
%%% fn is a text string specifying the location and name of the file to be
%%% written.
%%%
%%% type is the writing type: 'a' is append, 'w' is overwrite.
%%%
%%% sep is a string specifying the separater between fields in the output
%%% file (e.g.  ',' for csv files and '\t' for tab delimited files.
%%%
%%% Written by Aaron Schultz, May 5th, 2010 (aschultz@martinos.org)

if exist(fn)==0  || type == 'w'
    opt = 1;
else
    opt = 2;
end


if nargin < 4
    sep = '\t';
end

if nargin < 5
    if strcmpi(sep,',');
        pad = 0;
    else
        pad = 25;
    end
end

if nargin < 6
    percision = 5;
end

fid = fopen(fn, type);

if isstruct(S)
    if opt == 1
        flds = fields(S);
        for ii = 1:length(flds);
            if ii<length(flds);
                if pad > 0
                    fprintf(fid, ['%-' num2str(pad) 's'], flds{ii});
                else
                    fprintf(fid, '%s', flds{ii});
                end
            else
                if pad > 0
                    fprintf(fid, ['%-' num2str(pad) 's'], flds{ii});
                else
                    fprintf(fid, '%s', flds{ii});
                end
            end
            fprintf(fid, sep);
        end
        fprintf(fid, '\n');
    end
    
    for ii = 1:length(S);
        flds = fields(S(ii));
        for jj = 1:length(flds);
            ch = ischar(S(ii).(flds{jj})) || iscell(S(ii).(flds{jj}));
            if ch == 1
                if jj<length(flds)
                    if pad > 0
                        fprintf(fid, ['%-' num2str(pad) 's'], S(ii).(flds{jj}));
                    else
                        fprintf(fid, '%s', S(ii).(flds{jj}));
                    end
                else
                    if pad > 0
                        fprintf(fid, ['%-' num2str(pad) 's'], S(ii).(flds{jj}));
                    else
                        fprintf(fid, '%s', S(ii).(flds{jj}));
                    end
                end
            else
                if jj<length(flds)
                    if pad > 0
                        fprintf(fid,  ['%-' num2str(pad) '.' num2str(percision) 'f'], S(ii).(flds{jj}));
                    else
                        fprintf(fid,  ['%.' num2str(percision) 'f'], S(ii).(flds{jj}));
                    end
                else
                    if pad > 0
                        fprintf(fid,  ['%-' num2str(pad) '.' num2str(percision) 'f'], S(ii).(flds{jj}));
                    else
                        fprintf(fid, ['%.' num2str(percision) 'f'], S(ii).(flds{jj}));
                    end
                end
            end
            if jj<length(flds); fprintf(fid, sep); end
        end
        fprintf(fid, '\n');
    end
end

if iscell(S)
    if numel(size(S))>2
       error('CellArray Must be Two-Dimensional'); 
    end
    for ii = 1:size(S,1);
        for jj = 1:size(S,2);
            ch = ischar(S{ii,jj}) || iscell(S{ii,jj});
            if ch == 1
                if pad > 0
                    fprintf(fid, ['%-' num2str(pad) 's'], S{ii,jj});
                else
                    fprintf(fid, '%s', S{ii,jj});
                end
            else
                if pad > 0
                    fprintf(fid,  ['%-' num2str(pad) '.' num2str(percision) 'f'], S{ii,jj});
                else
                    fprintf(fid,  ['%.' num2str(percision) 'f'], S{ii,jj});
                end
            end
            if jj<size(S,2); fprintf(fid, sep); end
        end
        fprintf(fid, '\n');
    end
end

if isnumeric(S)
    if numel(size(S))>2
       error('Matrix Must be Two-Dimensional'); 
    end
    for ii = 1:size(S,1);
        for jj = 1:size(S,2);
            ch = ischar(S(ii,jj)) || iscell(S(ii,jj));
            if ch == 1
                if pad>0
                    fprintf(fid, ['%-' num2str(pad) 's'], S(ii,jj));
                else
                    fprintf(fid, '%s', S(ii,jj));
                end
            else
                if pad>0
                    fprintf(fid,  ['%-' num2str(pad) '.' num2str(percision) 'f'], S(ii,jj));
                else
                    fprintf(fid,  ['%.' num2str(percision) 'f'], S(ii,jj));
                end
            end
            if jj<size(S,2); fprintf(fid, sep); end
        end
        fprintf(fid, '\n');
    end
end


fclose(fid);