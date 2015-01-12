% Takes a number of brain images and runs an ANCOVA their voxelwise values
% with a particular set of X values (the parameter) and groups
% 
% 
function [groupF, groupP, parF, parP, interF, interP] = img2aoc(fileName, par, groups)

if ischar(fileName)
    fileName = strvcat2cell(fileName);
end

if length(par) ~= length(fileName)
    error('The number of images and the parameter are not equal...')
end

for f = 1:length(fileName)
    
    % Get image or matrix...
    if ischar(fileName{f})
        Y = spm_read_vols(spm_vol(fileName{f}));
    else
        Y = fileName{f};
    end
    
    % Linearise
    Y = Y(:);
    
    if f == 1
        V = spm_vol(fileName{f});
        data = nan(size(Y), length(par));         
        groupP = nan(V.dim);
        groupF = nan(V.dim);
        parP = nan(V.dim);
        parF = nan(V.dim);
        interP = nan(V.dim);
        interF = nan(V.dim);
    end

    data(:,f) = Y;
end

for v = 1:size(data,1)
    if all(and(isfinite(data(v,:)), data(v,:)~=0))
        
        [h,atab,ctab,stats] = aoctool(par, data(v,:), groups, 0.05, '', '', '', 'off');
        
        groupF(v) = atab{2,5};
        groupP(v) = atab{2,6};
        parF(v) = atab{3,5};
        parP(v) = atab{3,6};
        interF(v) = atab{4,5};
        interP(v) = atab{4,6};
    end    
end