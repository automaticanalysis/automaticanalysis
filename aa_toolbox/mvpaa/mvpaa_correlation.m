% MVPAA_CORRELATION - Simillate the betas/spmTs
% R - betas/spmTs (or residuals of betas/spmTs)

function [Simil] = mvpaa_correlation(aap, Resid)

% Rename settings to keep easier track...
EP = aap.tasklist.currenttask.settings;

%% R ==> (voxels, EP.conditions, EP.blocks, EP.sessions)]
sResid = reshape(Resid, [size(Resid,1), ...
    EP.conditions ...
    * EP.blocks ...
    * EP.sessions]);

% Set missing data to NaN here...
missed = all(sResid == 0, 1);
sResid(:,missed) = NaN;

% Simillate across voxels to find the correlation of voxel patterns
% across conditions.
if strcmp('Pearson', EP.corrType)
    % This is *much* faster than corr...
    Simil = corrcoef(sResid);
elseif strcmp('Spearman', EP.corrType)
    % Get Spearman correlations
    Simil = corr(sResid, 'type', 'Spearman');
elseif strcmp('Euclid', EP.corrType);
    % Get Euclidian distance
    Simil = squareform(pdist_complex(sResid', 'euclidean'));
elseif strcmp('sEuclid', EP.corrType);
    % Get Euclidian distance (standardised)
    Simil = squareform(pdist_complex(sResid', 'seuclidean'));
elseif strcmp('Mahalanobis', EP.corrType);
    % Get Mahalanobis distance
    dbstop if warning % If matrix is close to singular or badly scaled, we may see NaNs...
    Simil = squareform(pdist_complex(sResid', 'mahalanobis'));
else
    error('Incorrect metric of (dis)similarity between patterns');
end

if EP.triangulation ~= 2
    Simil = permute( ...
        reshape(Simil, ...
        [EP.conditions, ...
        EP.blocks * EP.sessions, ...
        EP.conditions, ...
        EP.blocks * EP.sessions]), ...
        [2, 4, 1, 3]);
end

%% Not that useful when betas are "solid". Default ignores it.
if EP.norm2z
    
    % Normalise the similarity measures within each block
    Simil = reshape(Simil, ...
        [(EP.blocks * EP.sessions)^2, ...
        (EP.conditions)^2]);
    Simil = (Simil - repmat(nanmean(Simil,2),[1, size(Simil,2)])) ...
        ./ repmat(nanstd(Simil,1,2), [1,size(Simil,2)]);
    Simil = reshape(Simil, ...
        [EP.blocks * EP.sessions, ...
        EP.blocks * EP.sessions, ...
        EP.conditions, EP.conditions]);
    
end