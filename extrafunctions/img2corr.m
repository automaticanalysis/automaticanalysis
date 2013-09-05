% Takes a number of brain images and correlates their voxelwise values with
% a particular set of values (the parmater)
function [C, T, P] = img2corr(fileName, par)

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
        C = nan(V.dim);
        P = nan(V.dim);
        T = nan(V.dim);
    end

    data(:,f) = Y;
end

for v = 1:size(data,1)
    if all(and(isfinite(data(v,:)), data(v,:)~=0))
        temp = corrcoef(data(v,:), par(:));
        C(v) = temp(2);
        [P(v), T(v)] = corr2pt(temp(2), length(par));
    end    
end