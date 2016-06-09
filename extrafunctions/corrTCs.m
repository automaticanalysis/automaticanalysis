%% This function correlates timecourses and displays a useful figure...
% Inputs:
%   aggrVars = PxC matrix, where P are timepoints across which to
%       correlate, and C are columns whose correlations with one another
%       you wish to see
%   TCnames = names for the columns
%   NaNnonsig = whether to NaN non-significant correlations (1 by default)
%   plotTCs = whether to plot the timecourses (1 by default)
function [sharedVar, h] = corrTCs(aggrVars, TCnames, NaNnonsig, plotTCs)

if nargin < 2
    aas_log([],true,'Not enough input variables')
end
if nargin < 3
    NaNnonsig = 1;
end
if nargin < 4
    plotTCs = 1;
end

%% SOME DEBUGGING CODE TO SHOW MOVE REGRESSOR INFO...

% And correlate the various regressors, as diagnostic
corrTC = corrcoef(aggrVars, 'rows', 'pairwise');

trigon = ~logical(tril(ones(size(corrTC))));
corrTC(~trigon) = NaN;

% Get the T and P values of the correlations...
pMR = nan(size(corrTC));
tMR = nan(size(corrTC));
dfV = nan(1,length(corrTC));

% Degrees of freedom for T and P values
for a = 1:length(corrTC);
    dfV(a) = sum(~isnan(aggrVars(:,a)));
end
dfMR = min(repmat(dfV', [1 length(dfV)]), repmat(dfV, [length(dfV) 1]));

[pMR tMR] = corr2pt(corrTC, dfMR);


% Bonferroni correct for multiple comparisons
pMR = pMR * sum(trigon(:));
% Get the correlations that are significant
ScorrTC = corrTC;
if NaNnonsig
    ScorrTC(pMR>0.05) = NaN;
end
sharedVar = sign(ScorrTC).*ScorrTC.^2;

MsharedVar = nanmean(abs(sharedVar(:)));

% Plot all the correlations that are significant on a figure

h = figure;
set(h, 'Position', [0 0 1200 700])

if plotTCs
    % Normalise the timecourses for display
    for t = 1:size(aggrVars,2)
        % Zero based
        aggrVars(:,t) = aggrVars(:,t) - min(aggrVars(:,t));
        % Max value is 1
        aggrVars(:,t) = aggrVars(:,t) ./ max(aggrVars(:,t));
    end
    
    subplot(1,2,1)
    %imagesc(aggrVars, 'AlphaData',~isnan(aggrVars))
    imagescnan(aggrVars)
    set(gca, 'Xtick', 1:length(corrTC), 'Ytick', [1 size(aggrVars,1)], ...
        'Xticklabel', TCnames, 'Yticklabel', [1 size(aggrVars,1)])
    ylabel('Timepoints')
    colorbar
    title('Variables across ''Timepoints''')
    rotateticklabel(gca,90);
    
    subplot(1,2,2)
end
%imagesc(sharedVar, 'AlphaData', ~isnan(sharedVar)) 
imagescnan(sharedVar)
set(gca, 'Xtick', 1:length(corrTC), 'Ytick', 1:length(corrTC), ...
    'Xticklabel', TCnames, 'Yticklabel', TCnames)
caxis([-1 1])
colorbar
axis equal tight
title(sprintf('Variance shared by variables. Mean: %0.2f %', MsharedVar*100))
rotateticklabel(gca,90);

pause(0.1)
end