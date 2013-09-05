% Takes a number of brain images and runs a linear model on their voxelwise values
% with a particular set of X values (the parameter) and groups
% 
% 
function [T, P, coefficientNames] = img2linearmodel(fileName, par, modelspec)

if nargin < 3
    modelspec = 'linear';
end

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
        
    end

    data(:,f) = Y;
end

testNum = [];
for v = 1:size(data,1)
    if all(and(isfinite(data(v,:)), data(v,:)~=0))
        
        mdl = LinearModel.fit(par, data(v,:), modelspec);
        tbl = anova(mdl);
        
        % We don't a priori know the size of our model, so...
        if isempty(testNum)
            testNum = length(mdl.CoefficientNames);
            coefficientNames = mdl.CoefficientNames;
            
            
            T = cell(testNum, 1);
            P = cell(testNum, 1);
            for p = 1:testNum
                T{p} = nan(V.dim);
                P{p} = nan(V.dim);
            end
        end
        
        for p = 1:testNum
            T{p}(v) = mdl.Coefficients.tStat(p);
            P{p}(v) = mdl.Coefficients.pValue(p);
        end
    end    
end