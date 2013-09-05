function cellR = uniqueCells(cellR)

% Get rid of repeated cell contents
for t1 = length(cellR):-1:1
    for t2 = t1-1:-1:1        
        if ischar(cellR{t1}) && ischar(cellR{t2})
            % If both cells are strings
            if strcmp(cellR{t1}, cellR{t2})
                cellR(t1) = [];
                break
            end            
        elseif isnumeric(cellR{t1}) && isnumeric(cellR{t2})
            % If both cells are numbers
            if cellR{t1}==cellR{t2}
                cellR(t1) = [];
                break
            end
        else
            % This function cannot cope with structures...
        end
    end
end