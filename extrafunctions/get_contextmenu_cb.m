function [h f] = get_contextmenu_cb(h,str)
h = get(h,'uicontextmenu');
Labels = textscan(str,'%s','delimiter','|'); Labels = Labels{1};

for i = 1:numel(Labels)
    ch = get(h,'children');
    l = get(ch,'label');
    h = ch(cell_index(l,Labels{i}));
end

f = get(h,'callback');