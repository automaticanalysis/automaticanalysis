% This function correlates timecourses and displays a useful figure...

function sharedVar = corrTCs(TCvals, TCnames)

%% SOME DEBUGGING CODE TO SHOW MOVE REGRESSOR INFO...

% Normalise the timecourses for display
for t = 1:size(TCvals,2)
    % Zero based
    TCvals(:,t) = TCvals(:,t) - min(TCvals(:,t));
    % Max value is 1
    TCvals(:,t) = TCvals(:,t) ./ max(TCvals(:,t));
end

% And correlate the various regressors, as diagnostic
corrTC = corrcoef(TCvals);

trigon = ~logical(tril(ones(size(corrTC))));
corrTC(~trigon) = NaN;
% Get the p and t values of the correlations...
[pMR tMR] = corr2pt(corrTC, size(TCvals,1));
% Bonferroni correct for multiple comparisons
pMR = pMR * sum(trigon(:));
% Get the correlations that are significant
ScorrTC = corrTC;
ScorrTC(pMR>0.05) = NaN;
sharedVar = sign(ScorrTC).*ScorrTC.^2;

MsharedVar = nanmean(abs(sharedVar(:)));

% Plot all the correlations that are significant on a figure

try close(2); catch; end
figure(2)
set(2, 'Position', [0 0 1200 700])

subplot(1,2,1)
imagesc(TCvals)
set(gca, 'Xtick', 1:length(corrTC), 'Ytick', [1 size(TCvals,1)], ...
    'Xticklabel', TCnames, 'Yticklabel', [1 size(TCvals,1)])
ylabel('Timepoints')
colorbar
title('Variables across ''Timepoints''')
rotateticklabel(gca,90);

subplot(1,2,2)
imagescnan(sharedVar)
set(gca, 'Xtick', 1:length(corrTC), 'Ytick', 1:length(corrTC), ...
    'Xticklabel', TCnames, 'Yticklabel', TCnames)
caxis([-1 1])
colorbar
axis equal
title(sprintf('Variance shared by variables. Mean: %0.2f %', MsharedVar*100))
rotateticklabel(gca,90);

pause(0.1)
end