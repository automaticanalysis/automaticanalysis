function [ostr, rem] = strtok_ptrn(istr,ptrn)
ostr = istr; rem = '';
isfx = strfind(istr,ptrn);
if ~isempty(isfx)
    ostr = istr(1:isfx(1)-1);
    rem = istr(isfx(1):end);
end