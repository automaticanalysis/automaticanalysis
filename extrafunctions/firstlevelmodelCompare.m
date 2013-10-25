%% This function compares the output of SPM models at the 1st level
% Inputs:
%   pths = paths to your models in a cell array
%   maskVolVol = which mask to use to mask in your voxels of interest
%       (defaults to using all voxels that are finite and not 0);

function firstlevelmodelCompare(pths, maskVol)
if nargin < 1
    error('Not enough inputs, need paths to both model A and model B')
elseif length(pths) == 1
    error('Not enough models for comparison')
end

if nargin < 2
    maskVol = [];
end

% Get maskVol?
if ~isempty(maskVol)
    if ischar(maskVol)
        maskVol = spm_read_vols(spm_vol(maskVol));
    end
    maskVol = logical(maskVol);
end
    
for p = 1:length(pths);
    % Get betas, SPM for model A
    D{p} = dir(fullfile(pths{p}, 'beta*img'));
    load(fullfile(pths{p}, 'SPM.mat'));
    SPMcols{p} = SPM.xX.iC;
    SPMnames{p} = SPM.xX.name;
        
    % INTEREST AND NUISANCE COLUMNS...
    SPMinterest{p} = SPM.xX.name(SPM.xX.iC);
    SPMnuisance{p} = SPM.xX.name(SPM.xX.iG);

    % Check...
    if p > 1
        if length(SPMinterest{p}) ~= length(SPMinterest{p-1})
            error('Your models have a different number of interest columns')
        end
        for d = 1:length(SPMinterest{p})
           if ~strcmp(SPMinterest{p}{d}, SPMinterest{p-1}{d})
               warning('Your models have a condition names')
           end
        end
    end
end

%% CORRELATIONS BETWEEN BETAS...
h.Fig = figure;
splitFig = ceil(sqrt(length(SPMinterest{1})));

C = cell(size(SPMinterest{1}));

for d = 1:length(SPMinterest{1})
    for p = 1:length(pths)
        V = spm_vol(fullfile(pths{p}, D{p}(SPMcols{p}(d)).name));
        Y = spm_read_vols(V);
        % maskVol things we don't want...
        if ~isempty(maskVol);
            Y = Y.*maskVol;
        end
        Y = Y(isfinite(Y) & Y~=0);
        if p == 1
            data = nan(length(Y), length(pths));
        end
        data(:,p) = Y;
    end
    
    % Make correlation matrix for models
    C{d} = corrcoef(data, 'rows', 'pairwise');
    C{d} = triu(C{d},1);
    C{d}(C{d}==0) = NaN;
    
    % Plot
    subplot(splitFig, splitFig, d)
    imagescnan(C{d});
    caxis([-1, 1]);
    axis square off
    title(SPMinterest{1}{d})
end
