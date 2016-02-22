function out_str = list_index(line,count,subcount,shave)
if nargin < 4, shave = 1; end % shave tailing _
if nargin < 3, subcount = 1; end
out_str = '';
for s = 1:numel(subcount)
    val_str = sscanfitem(strrep(sscanfitem(strrep(line,';',' '),count),'_',' '),subcount(s));
    % sublist = list_item(line,count,';'); % BASH-like behavior
    % val_str = list_item(sublist,subcount,'_');
    if ~isempty(val_str), out_str = [out_str val_str '_']; end
end
if ~isempty(out_str) && shave, out_str = out_str(1:end-1); end
end

% function str = list_item(line,ind,sep)
% line = strrep(line,sep,' ');
% str = sscanfitem(line,ind);
% while isempty(str)
%     ind = ind - 1;
%     str = sscanfitem(line,ind);
% end
% end