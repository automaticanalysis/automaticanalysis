% aamod_clusteringXValidation
%
% Details to come....

% CW - 2014-03-17

function [aap, resp] = aamod_clusteringXValidation(aap, task)

resp='';

switch task
    case 'report'
        
    case 'doit'
        
        numSubjs = length(aap.acq_details.subjects);
        
        % Get the settings for this task
        settings = aap.tasklist.currenttask.settings;
        thisModNum = aap.tasklist.currenttask.modulenumber;
        
        % We need info about the inputs... for double checking ROIs, etc.
        inputSources = aap.internal.inputstreamsources{thisModNum};
        srcModName = inputSources.stream(1).sourcestagename;
        srcModNum = inputSources.stream(1).sourcenumber;
        
        % Get the ROI from task settings; it could be in two places :P
        try
            whichROI = aap.tasklist.main.module(srcModNum).extraparameters.aap.tasklist.currenttask.settings.whichROI;
        catch
            srcModInd = aap.tasklist.main.module(srcModNum).index;
            srcModSettings = aap.tasksettings.(srcModName)(srcModInd);
            whichROI = srcModSettings.whichROI;
        end
        
        % Empty struct array for initialization...
        connectivityMatrices = struct('seed', {}, 'seedVoxMm', {}, 'seedVoxInd', {}, 'seedSpace', {}, 'targetNames', {}, 'targetVoxInd', {}, 'targetSpace', {}, 'correlationMatrix', {});

        % Collect the connectivity matrices for each subject
        for subInd = 1:numSubjs
            matrixFile = aas_getfiles_bystream(aap, subInd, 'firstlevel_fConnMatrix_avg');
            load(matrixFile);
            
            if ~exist('avgMatrices', 'var');
                aas_log(aap, 1, sprintf('%s doesn''t contain a data structure called ''avgMatrices''', matrixFile));
            end
            
            % Could crash here if a subject has a different number of ROIs
            % than the others, but I'm too lazy to check, so instead you
            % can just read this line if it crashes here.
            connectivityMatrices(end+1, :) = avgMatrices; 
        end
        
        numROIs = size(connectivityMatrices, 2);
        
        % Which ROI do we want analyze?
        if isnumeric(whichROI)            
            roiI = whichROI;    % Could be set using an index...
            if roiI > numROIs, aas_log(aap, 1, 'There are only %d seed ROIs, but whichROI=%d', numROIs, roiI); end 
        
        else                    % Could be specified as a string...
            roiI = -1;
            for i = 1 : numROIs
               roiNames = unique({connectivityMatrices(:,i).seed});              
               m = regexp(roiNames{1}, whichROI);
               if m
                   roiI = i;
                   break
               end 
            end
            if roiI < 0, aas_log(aap, 1, sprintf('Can''t find ROI with name %s', whichROI)); end
        end
        
        % Check the ROI names, warn user if they don't correspond across subjects
        roiNames = unique({connectivityMatrices(:,roiI).seed});      
        if length(unique(roiNames)) > 1
            aas_log(aap, 0, sprintf('\nWarning: seed ROIs might be different across subjects.\nUsing these seeds: %s\n', strjoin(roiNames, ', ')), 'red');
        else
            aas_log(aap, 0, sprintf('\nWorking with seed region: %s', roiNames{1}));
        end       
        
        % Trim the ROIs we aren't using
        connectivityMatrices = connectivityMatrices(:, roiI);
        
        
        % Let's get the group parcellation, so we can find correspondance
        % among the individuals
        groupInfo = load(aas_getfiles_bystream(aap, 'group_module_info'));
                
        % This is the number of modules we are going to deal with
        numModules = sum([groupInfo.moduleInfo.numVox] > settings.minVoxPerMod);
        
        % Other important numbers!
        numSeedVox = size(groupInfo.moduleInfo(1).spatialPattern, 2);
        numTargets = size(groupInfo.moduleInfo(1).connectivityPattern, 2);
        
        % Get the spatial patterns of the top N modules
        groupNSpatialPatterns = vertcat(groupInfo.moduleInfo(1:numModules).spatialPattern);
        aas_log(aap, 0, sprintf('Cross-Validating spatial and connectivity profiles for %d modules, across %d subjects', numModules, numSubjs));
       
        % Now build the spatial connectivity fingerprints for each subject using their LOO maps
        subjConnectivityPatterns = zeros(numModules, numTargets, numSubjs);
        subjSpatialPatterns = zeros(numModules, numSeedVox, numSubjs);
        
        for subjInd = 1 : numSubjs
            
            fprintf('Working with LOO parcellations generated for %s\n', aap.acq_details.subjects(subjInd).mriname);
            
            % Get the LOO generated modules
            looInfo = load(aas_getfiles_bystream(aap, subjInd, 'N-1_module_info'));
            looSpatialPatterns = vertcat(looInfo.moduleInfo.spatialPattern);
            
            % Spatial similiarity between top N group modules, and all of
            % this individual's LOO modules
            spatialSimil = corr(groupNSpatialPatterns', looSpatialPatterns');
            
            % Find the mappings, which group module maps to which loo
            % module ('loo' refers to the data from leave-one-out, 'subj' 
            % reconstructed individual data). This is hard because because
            % we need a unique one-to-one mapping
            group2loo = []; c = [];
            for m = 1 : numModules
                
                % Module pair with greatest similarity
                [c(m), i] = max(spatialSimil(:));
                
                % These indices are the mapping
                [group2loo(m,1), group2loo(m,2)] = ind2sub(size(spatialSimil), i);
                
                % Now remove those modules from the similarity matrix
                spatialSimil(group2loo(m,1),:) = -Inf;
                spatialSimil(:,group2loo(m,2)) = -Inf;
                
            end % end for each module
            
            % It's not uncommon for some modules to have NO overlap?  What
            % do we do??
            if any(c<=0)
                aas_log(aap, 0, sprintf('Warning: %d of the group modules don''t overlap with any in this LOO iteration! Try increasing voxel threshold??', sum(c<=0)));
            end
            
            % Sort the mappings so they make sense again
            [junk,i] = sort(group2loo(:,1), 'ascend');
            group2loo = group2loo(i,:);
            
            % Now let's reconstruct the connectivity fingerprints of each module in the left out individual
            for modI = 1 : numModules
                
                % Which module index in the loo data?
                looModI = group2loo(modI, 2);
                
                % Voxel indices for this loo module 
                modVoxI = find(looInfo.moduleInfo(looModI).spatialPattern);
                
                % Reconstruct the connectivity profile of this module by
                % averaging the profile across all voxels in this subject 
                % that belong to the LOO module.
                subjConnectivityPatterns(modI, :, subjInd) = nanmean(connectivityMatrices(subjInd).correlationMatrix(modVoxI, :));                              
                
            end
            
            % Next is harder, let's reconstruct the spatial pattern for
            % each module by assigning voxels to modules based on their
            % connectivity profiles.    
            moduleConnectivity = vertcat(looInfo.moduleInfo(group2loo(:,2)).connectivityPattern)'; % Clever, reshuffle for group modulemembership here!
            seedVoxConnectivity = connectivityMatrices(subjInd).correlationMatrix';
            
            % Filter NaNs
            keepTargets = ~any(isnan([seedVoxConnectivity moduleConnectivity]),2);
            
            % Compute the similiarity between each voxel and each module
            connectivitySimil = corr(seedVoxConnectivity(keepTargets,:), moduleConnectivity(keepTargets,:), 'type', 'Spearman');
            
            % Find the most similar module for each voxel
            [junk, bestMod] = max(connectivitySimil, [], 2);
            
            % Create the binary spatial pattern, save for this subject
            subjSpatialPatterns(:,:,subjInd) = (repmat(bestMod, 1, numModules) == repmat(1:numModules, numSeedVox, 1))';
            
        end % end for each subject
        
        % Create and save probability maps for each module
        modulePath = aas_getpath_bydomain(aap, aap.tasklist.currenttask.domain, []);
        imagePath = fullfile(modulePath, 'pmaps');
        if ~exist(imagePath), mkdir(imagePath); end
        
        space = connectivityMatrices(1).seedSpace;
        
        seedVoxI = connectivityMatrices(1).seedVoxInd;
        
        % Create SPM Vol structure to use as a template for our module maps
        Vt = struct('fname', '', ...
                    'dim',   space.dim, ...
                    'mat',   space.mat, ...
                    'dt',    [4 0], ...      % Don't need much storage for these images
                    'pinfo', [1 0 0]', ...
                    'n',     [1 1]);
        pMaps = [];

        fprintf('Writing out probability maps...');
        for m = 1 : numModules

            % Probability of voxel belonging to this module, across
            % subjects
            pmod = mean(subjSpatialPatterns(m,:,:), 3);
            
            % Scale to fit the data type 
            pmod = single(pmod * 2^16);
            
            Ym = zeros(space.dim);
            Ym(seedVoxI) = pmod;
            
            Vm = Vt;
            Vm.fname = fullfile(imagePath, sprintf('MOD%02d_pmap.nii', m));
            spm_write_vol(Vm, Ym);
            
            pMaps = char(pMaps, Vm.fname);
            
        end
        
        pMaps(1,:) = [];
        aap = aas_desc_outputs(aap, [], 'XVal_mod_pmaps',pMaps);   
        
        fprintf('done');
        
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap, 1, sprintf('Unknown task %s',task));
end

end
