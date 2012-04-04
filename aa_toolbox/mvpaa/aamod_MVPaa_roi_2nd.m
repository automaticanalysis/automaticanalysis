% AA module - MVPaa 2nd level
%
% Modified for aa4 by Alejandro Vicente-Grabovetsky Feb-2011

function [aap,resp] = aamod_MVPaa_roi_2nd(aap,task)

resp='';

switch task
    case 'doit'
        
        tic
        
        Stats = []; meanSimil = [];% Statistic structure that we load for each participant
        for p = 1:length(aap.acq_details.subjects)
            load(aas_getfiles_bystream(aap,p,'MVPaa'));
            
            if p == 1
                % MVPA data for Multivariate results
                indStats = nan(size(Stats,1), ...
                    length(aap.acq_details.subjects), ...
                    size(Stats,2), ...
                    size(Stats,3));
                
                % MVPA data in matrix form..
                aggrSimil = nan(size(meanSimil,1), ...
                    length(aap.acq_details.subjects), ...
                    size(meanSimil,2), ...
                    size(meanSimil,3));
                
                % Set some of the settings in 2nd level from 1st level...
                aap.tasklist.currenttask.settings.tests = EP.tests;
            end
            
            % Gather data from each participant
            indStats(:, p, :,:) = Stats;
            aggrSimil(:, p, :,:) = meanSimil;
        end
        
        % 2nd level stats (aggregate StatisticollStat -> collapsed Statistics)
        aggrStats = mvpaa_statistics_2nd(aap, indStats);
        
        % Mean of the meanSimil data
        aggrSimil = squeeze(nanmean(aggrSimil,2));
        
        %% Describe outputs...
        save(fullfile(aap.acq_details.root, 'MVPaa.mat'), ...
                'indStats', 'aggrStats', 'aggrSimil', 'EP')
        aap=aas_desc_outputs(aap,'MVPaa_2nd', fullfile(aap.acq_details.root, 'MVPaa.mat'));
        
        time_elapsed
end