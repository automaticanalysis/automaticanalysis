function [in ord] = matchAll(pl, rv, trm);

if nargin == 2;
    trm = [0 0];
end

in = [];
ord = [];

c = 0;
for ii = 1:length(pl);
    tmp = strmatch(pl{ii}(1+trm(1):end-trm(2)),rv);
%     tmp = contains(pl{ii}(1+trm(1):end-trm(2)),rv);
    if ~isempty(tmp)
        c = c+1;
        in(c) = ii; 
%         try
            ord(c) = tmp;
%         catch
%             keyboard;
%         end
    end
end