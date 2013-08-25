% MVPAa_STATISTINGS_S Sort how to do statistics
% SU = Subject statistics for each ROI, contrast and test
% contrasts = contrasts used

function aggrStats = mvpaa_statistics_2nd(aap, indStats)

%% Statistics...
if strcmp(aap.tasklist.currenttask.settings.statsType2nd, 'ttest')
    aggrStats = nan(size(indStats,1), size(indStats,3), 5);
elseif strcmp(aap.tasklist.currenttask.settings.statsType2nd, 'signrank')
    aggrStats = nan(size(indStats,1), size(indStats,3), 3);
else
    aas_log(aap,true, sprintf('No such 2nd level statistic test as: %s', ...
        aap.tasklist.currenttask.settings.statsType2nd));
end

%SU => ( EP.rois, EP.participants, size(contrasts,2), length(EP.tests));
indStats = indStats(:,:,:,1); % Either beta, mean or median; don't care about t-value or p-value...

for r = 1:size(indStats,1) % ROIs
    for c = 1:size(indStats,3) % Contrasts
        % Remove nans from data...
        data = squeeze(indStats(r,:,c));
        data(isnan(data)) = [];
        
        if strcmp(aap.tasklist.currenttask.settings.statsType2nd, 'ttest')
            aggrStats(r,c,1) = mean(data);
            [junk, p, junk, stats] = ttest(data);
            aggrStats(r,c,2) = stats.tstat;
            aggrStats(r,c,3) = p;
            aggrStats(r,c,4) = stats.sd / sqrt(stats.df);
            % JB test on all values of correlations? (or absolute values?)
            warning off
            [junk, aggrStats(r,c,5)] = jbtest(data);
            warning on
        elseif strcmp(aap.tasklist.currenttask.settings.statsType2nd, 'signrank')
            aggrStats(r,c,1) = median(data);
            % Signed Rank for one samples
            p = signrank(data);
            aggrStats(r,c,2) = tinv(1-p, length(data)-1);
            aggrStats(r,c,3) = p;
        end
    end
end