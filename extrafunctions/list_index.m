function val_str = list_index(line,count,subcount)
if nargin < 3, subcount = 1; end
val_str = sscanfitem(strrep(sscanfitem(strrep(line,';',' '),count),'_',' '),subcount);
% sublist = list_item(line,count,';'); % BASH-like behavior
% val_str = list_item(sublist,subcount,'_');
end

function str = list_item(line,ind,sep)
line = strrep(line,sep,' ');
str = sscanfitem(line,ind);
while isempty(str)
    ind = ind - 1;
    str = sscanfitem(line,ind);
end
end