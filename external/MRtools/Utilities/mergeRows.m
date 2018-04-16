function td = mergeRows(dat,criteria)
%%
ind = [];
if ~isnumeric(criteria) && ~isobject(criteria)
    criteria = nominal(criteria);
end
list = unique(criteria,'stable');


cols = dat.Properties.VariableNames;

for ii = 1:numel(list)
    i1 = find(criteria==list(ii));
    
    for jj = 1:numel(cols)
        if isnumeric(dat.(cols{jj}))
            dat.(cols{jj})(i1) = nanmean(dat.(cols{jj})(i1));
        end
    end
    
    ind(ii,1) = i1(1);
end

td = dat(ind,:);
for jj = 1:numel(cols)
    if isobject(td.(cols{jj}))
        tmp = cellstr(td.(cols{jj}));
        td.(cols{jj}) = nominal(regexprep(tmp,'<undefined>',''));
    end
end



