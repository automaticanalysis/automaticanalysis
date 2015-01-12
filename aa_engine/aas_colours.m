function colours = aas_colours(colnum)
if nargin < 1
    colnum = 15;
end

% R G B Y V C
% rK, rG, rB,
try    
    colours = distinguishable_colors(colnum);
    colours = mat2cell(colours, ones(1,colnum), 3);
catch
    colours = {[1 0 0], [0 1 0], [0 0 1], [0 1 1], [1 0 1], [1 1 0], ...
        [1 0.5 0.5], [0.5 1 0.5], [0.5 0.5 1], ...
        [1 0.5 0], [1 0 0.5], [0.5 1 0], [0.5 0 1], [0 1, 0.5], [0 0.5 1]};
end

