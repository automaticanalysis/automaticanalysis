% AA module - MVPaa 1st level (ROI based)
%
% Modified for aa4 by Alejandro Vicente-Grabovetsky Feb-2011

function [aap,resp] = aamod_MVPaa_roi_1st(aap,task,p)

resp='';

switch task
    case 'doit'
        
        %% PREPARATIONS
        
        subjname = strtok(aap.acq_details.subjects(p).subjname, '/');
        fprintf('Working with data from participant %s. \n', subjname)
        
        % Get the contrasts for this subject...
        aap.tasklist.currenttask.settings.contrasts = mvpaa_loadContrasts(aap,p);
                       
        % Load the ROIs from which to extract the data
        ROIimg = aas_getfiles_bystream(aap,p,'rois');
        
        % Which tests will we use?
        if ~isempty(findstr(aap.tasklist.currenttask.settings.statsType, 'GLM'))
            aap.tasklist.currenttask.settings.tests = {'beta', 't-value', 'p-value', 'SE'};
        elseif ~isempty(findstr(aap.tasklist.currenttask.settings.statsType, 'ttest'))
            aap.tasklist.currenttask.settings.tests = {'mean', 't-value', 'p-value', 'SE', 'normality'};
        elseif ~isempty(findstr(aap.tasklist.currenttask.settings.statsType, 'signrank'))
            aap.tasklist.currenttask.settings.tests = {'median', 't-value (est)', 'p-value'};
        end        
        
        %% ANALYSIS
        
        % Load the data into a single big structure...
        [aap data] = mvpaa_loadData(aap, p);
        
        ROInum = size(ROIimg,1);
        
        % Create output arrays...
        Stats = NaN(ROInum, ...
            length(aap.tasklist.currenttask.settings.contrasts), ...
            length(aap.tasklist.currenttask.settings.tests));
        meanSimil = NaN(ROInum, ...
            aap.tasklist.currenttask.settings.conditions, ...
            aap.tasklist.currenttask.settings.conditions);
        
        % Loop the routine over all ROIs
        for r = 1:ROInum
            [Rpth Rfn Rext] = fileparts(deblank(ROIimg(r,:)));
            
            % Extract betas from each: ROI, voxel, condition, subblock, session
            ROI = int8(spm_read_vols(spm_vol(fullfile(Rpth, [Rfn Rext]))));
            
            % Check that the ROI size is equal to the data size
            if any(size(ROI) ~= size(data{1,1,1}));
                aas_log(aap, true, 'Your ROI size is different from your data size!');
            end
            
            % Trick for non-binary ROIs...
            if length(unique(ROI))>2
                ROI = ROI > 0;
            end
            voxels = sum(ROI(:));
            
            % ROI to linear index...
            ROI = find(ROI);
            
            Betas = mvpaa_extraction(aap, data, ROI, voxels);
            
            fprintf('\t ROI = %s; vox. = %d (%d)\n',Rfn, sum(~isnan(data{1,1,1}(ROI))), voxels)
            
            if isempty(Betas)
                aas_log(aap, false, sprintf('Not enough voxels in ROI, minimum is %i, you have %i', ...
                    aap.tasklist.currenttask.settings.minVoxels, sum(~isnan(data{1,1,1}(ROI)))));
                continue
            end
            
            % Get the residuals
            Resid = mvpaa_shrinkage(aap, Betas);
            
            % Compute similarities of the the data
            Simil = mvpaa_correlation(aap, Resid);
            
            [Stats(r,:,:), meanSimil(r, :,:)] = mvpaa_statistics(aap, Simil);
            
            %% SOME DIAGNOSTICS...
            if aap.tasklist.currenttask.settings.diagnostic
                if ~exist(fullfile(aap.acq_details.root, 'diagnostics'), 'dir')
                    mkdir(fullfile(aap.acq_details.root, 'diagnostics'))
                end
                
                try close(2); catch;end
                figure(2)
                set(2, 'Position', [0 0 1200 500], 'Name', Rfn)
                
                subplot(1,3,1)
                imagesc( ...
                    reshape(Resid, ...
                    [size(Resid,1) size(Resid,2)*size(Resid,3)*size(Resid,4)]))
                axis equal
                axis off
                title('Residuals')
                
                subplot(1,3,2)
                imagesc( ...
                    reshape(permute( ...
                    Simil, [3, 1, 4, 2]), ...
                    [size(Simil,3)*size(Simil,1), size(Simil,4)*size(Simil,2)]));
                caxis([-1 1])
                axis equal
                axis off
                title('Similarity matrix...')
                
                subplot(1,3,3)
                imagesc(squeeze(meanSimil(r, :,:)));
                caxis([-1 1])
                axis equal
                axis off
                title('...collapsed across sessions and blocks')
                
                print('-djpeg','-r75',fullfile(aap.acq_details.root, 'diagnostics', ...
                [mfilename '__' subjname '_' num2str(r) '.jpeg']));
            end     
            try close(2); catch;end
        end
        
        %% DESCRIBE OUTPUTS
        EP = aap.tasklist.currenttask.settings;
        save(fullfile(aas_getsubjpath(aap,p), [subjname '.mat']), ...
            'meanSimil', 'Stats', 'EP')
        aap=aas_desc_outputs(aap,p,'MVPaa', fullfile(aas_getsubjpath(aap,p), [subjname '.mat']));
end