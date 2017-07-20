function [ostr, rem] = strtok_ptrn(istr,ptrn)
ostr = istr; rem = '';
isfx = strfind(istr,ptrn);
if iscell(isfx)
    for i = 1:numel(istr)
        if ~isempty(isfx{i})
            ostr{i} = istr{i}(1:isfx{i}(1)-1);
            rem{i} = istr{i}(isfx{i}(1):end);
        else
            rem{i} = '';
        end
    end
else
    if ~isempty(isfx)
        ostr = istr(1:isfx(1)-1);
        rem = istr(isfx(1):end);
    end    
end