%aamod_fconn_computematrix Compute functional connectivity matrix
%
%  aamod_fc_computematrix(aap, subjInd, task) computes a functional
%  connectivity matrix based on timecourse relationships across voxels.
function [aap, resp] = aamod_fconn_computematrix(aap, task, subjInd)

resp='';

switch task
    case 'report'
        
    case 'doit'                        
        settings = aap.tasklist.currenttask.settings;
        
        fconn = [];
        fconn.measure = 'pearson'; % to be used later?
        
        
        %% If voxelwise, read in residual images
        
        % Get image names
        spmName = aas_getfiles_bystream(aap, subjInd, 'firstlevel_spm');
        maskName = aas_getfiles_bystream(aap, subjInd, 'firstlevel_brainmask');
        roiName = aas_getfiles_bystream(aap, subjInd, 'roi'); % source voxels
        epiNamesAll = aas_getfiles_bystream(aap, subjInd, 'epi'); % typically filtered, or first level residuals
        
        
        % Just keep img/nii files (no hdr)
        epiNames = {};
        for i = 1:size(epiNamesAll,1)
            thisName = strtok(epiNamesAll(i,:));
            if ismember(thisName(end-3:end), {'.img' '.nii'})
                epiNames{end+1} = thisName;
            end
        end
        
        % read in SPM.mat
        load(spmName);
        
        % Make sure the number of images we found corresponds to those listed in
        % SPM.nscan
        if length(epiNames)~=sum(SPM.nscan(aap.acq_details.selected_sessions))
            aas_log(aap, true, sprintf('You don''t have the right number of timeseries images (expected %d, found %d).', sum(SPM.nscan(aap.acq_details.selected_sessions)), length(epiNames)));
        end
        
        
        % Set up target voxels (typically brain mask)
        [brainMask, XYZtarget] = spm_read_vols(spm_vol(maskName));
        targetInd = find(brainMask>0);
        nTargetVoxels = length(targetInd);
        XYZtarget = XYZtarget(targetInd);
        
        
        % Set up source voxels (only selected, or by default, same as
        % target for a symmetrical matrix)
        [sourceMask, XYZsource] = spm_read_vols(spm_vol(roiName));
        sourceInd = find(sourceMask>0);
        nSourceVoxels = length(sourceInd);
        XYZsource = XYZsource(sourceInd);
        
        
        
        aas_log(aap, false, sprintf('Calculating connectivity between %d source voxels and %d target voxels.', nSourceVoxels, nTargetVoxels));
                                
        
        
        % keep track of correlations and p values
        
        aas_log(aap, false, 'Making big maatrix to hold correlations...');        
        fconn.rSingle = zeros(nSourceVoxels, nTargetVoxels, 'single'); % NB single
        aas_log(aap, false, 'Matrix made.');
        
        
        % Do correlations across all volumes (concatenate) or split by
        % session and then average?
        if settings.concatenate
            V = spm_vol(strvcat(epiNames));            
            [Y, XYZ] = spm_read_vols(V);
            sizeY = size(Y);
            
            tic
            
            for s = 1:nSourceVoxels
                
                fprintf('Source voxel %d/%d...', s, nSourceVoxels);
                
                [si, sj, sk, sz] = ind2sub(size(sourceMask), sourceInd(s));
                sourceData = Y(si,sj,sk,:);
                
                % If we didn't estimate a source voxel, don't do it;
                % otherwise, do it.
                if ~isfinite(sourceData(1))
                    fconn.rSingle(s,:) = NaN;
                else
                    
                    for t = 1:nTargetVoxels
                        [ti, tj, tk, tz] = ind2sub(size(brainMask), targetInd(t));
                        
                        targetData = Y(ti,tj,tk,:);
                        
                        [r,p] = corrcoef(sourceData, targetData);
                        fconn.rSingle(s,t) = single(r(1,2));
                        
                    end % target
                end % checking for source NaN
                
                hElapsed = toc/60/60;
                fprintf('done. %.2f%% done in %.2f hours (estimate %.2f h remaining).\n', 100*s/nSourceVoxels, hElapsed, (hElapsed/s)*(nSourceVoxels-s));
                
            end % source            
        else            
            % For each session, read in data and do correlations. We'll sum
            % data in each matrix location to avoid storing multiple copies of
            % the matrix, and then divide by num sessions at the end.
            for thisSess = aap.acq_details.selected_sessions
                
                aas_log(aap, false, sprintf('Processing session %d...', thisSess));
                
                % Figure out which scans belong to this session
                if thisSess==1
                    sessIndex = 1:SPM.nscan(thisSess);
                else
                    sessIndex = sum(SPM.nscan(1:thisSess-1))+1:sum(SPM.nscan(1:thisSess));
                end
                                
                % read in images
                V = spm_vol(strvcat(epiNames(sessIndex)));
                [Y, XYZ] = spm_read_vols(V);
                sizeY = size(Y);                
                
                for s = 1:nSourceVoxels
                    [si, sj, sk, sz] = ind2sub(size(sourceMask), sourceInd(s));
                    sourceData = Y(si,sj,sk,:);
                    
                    % If we didn't estimate a source voxel, don't do it;
                    % otherwise, do it.
                    if ~isfinite(sourceData(1))
                        fconn.rSingle(s,:) = NaN;
                    else
                        
                        for t = 1:nTargetVoxels
                            [ti, tj, tk, tz] = ind2sub(size(brainMask), targetInd(t));
                            
                            targetData = Y(ti,tj,tk,:);
                            
                            [r,p] = corrcoef(sourceData, targetData);
                            fconn.rSingle(s,t) = single(double(fconn.rSingle(s,t)) + r(1,2));
                            
                        end % target
                    end % checking for source NaN
                end % source
                
            end % looping through sessions
            
            % average over sessions for r values
            nSessions = numel(aap.acq_details.selected_sessions);
            fconn.rSingle = fconn.rSingle ./ nSessions;
        end
            
        
        % output
        [pth, nm, ext] = fileparts(spmName);
        subName = aap.acq_details.subjects(subjInd).mriname;
        fconnName = sprintf('%s-fconn%s.mat', subName, settings.matrixsuffix);
        fconnPath = fullfile(pth, fconnName);
        
        % Add descriptives
        fconn.XYZsource = XYZsource;
        fconn.XYZtarget = XYZtarget;
        fconn.statistic = 'Pearson r';
        fconn.space = 'MNI';
        fconn.version = '00.0';
        fconn.path = fconnPath;
        fconn.README = 'First iteration of connectivity matrix. What''s the greek letter before alpha?';                
        
        save(fconnPath, 'fconn');
        
        % If there's a matrix farm, make a link
        if ~isempty(settings.matrixfarm)
            try
                if ~isdir(settings.matrixfarm)
                    mkdir(settings.matrixfarm);
                end
                
                farmPath = fullfile(settings.matrixfarm, fconnName);                
                
                % Copy or link, depending
                if settings.matrixfarmlink==0
                    system(sprintf('cp %s %s', fconnPath, farmPath));
                else
                    system(sprintf('ln -s %s %s', fconnPath, farmPath));
                end
            catch
                aas_log(aap, false, sprintf('Tried to link .mat file to matrixfarm %s but it did not work.', settings.matrixfarm));
            end
        end
            
        
        aap = aas_desc_outputs(aap, subjInd, 'firstlevel_fconn_matrix_mat', fconnPath);
        
        
        
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end




