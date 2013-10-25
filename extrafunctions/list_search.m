function ind = list_search(line,patt)
i = 1;
while ~isempty(list_index(line,i))
    if strcmp(list_index(line,i), patt), break; end
    i = i + 1;
end
ind = ~isempty(list_index(line,i)) * i;
