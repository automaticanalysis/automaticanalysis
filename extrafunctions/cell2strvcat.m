function vCatStr = cell2strvcat(strCell)
    if ~iscell(strCell)
        error('Input is not a cell array')
    end    
    vCatStr = '';
    for o = 1:length(strCell)
        vCatStr = strvcat(vCatStr, strCell{o});
    end
end