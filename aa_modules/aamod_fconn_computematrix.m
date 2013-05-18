%aamod_fconn_computematrix Compute functional connectivity matrix
%
%  aamod_fc_computematrix(aap, subjInd, task) computes a functional
%  connectivity matrix based on timecourse relationships across voxels.
function [aap, resp] = aamod_fconn_computematrix(aap, task, subjInd)

resp='';

switch task
    case 'report'
        
    case 'doit'        
        
        fconn = [];
        fconn.measure = 'pearson';
        
        
        %% If voxelwise, read in residual images
        
        % Get image names
        spmName = aas_getfiles_bystream(aap, subjInd, 'firstlevel_spm');
        maskName = aas_getfiles_bystream(aap, subjInd, 'firstlevel_brainmask');
        %betaNames = aas_getfiles_bystream(aap, subjInd, 'firstlevel_betas');
        residNamesAll = aas_getfiles_bystream(aap, subjInd, 'firstlevel_residuals');
        
        % Just keep img/nii files (no hdr)
        residNamesImg = {};
        for i = 1:size(residNamesAll,1)
            thisName = strtok(residNamesAll(i,:));
            if ismember(thisName(end-3:end), {'.img' '.nii'})
                residNamesImg{end+1} = residNamesAll(i,:);
            end
        end
        
        % read in SPM.mat
        load(spmName);
        
        % Make sure the number of images we found corresponds to those listed in
        % SPM.nscan
        if length(residNamesImg)~=sum(SPM.nscan)
            aas_log(aap, true, 'You don''t have the right number of residual images.');
        end
        
        % Read in the brain mask and only use voxels that were estimated
        [brainMask, XYZ] = spm_read_vols(spm_vol(maskName));
        goodVoxelInd = find(brainMask>0);
        nVoxels = length(goodVoxelInd);
        XYZ = XYZ(goodVoxelInd); % should be the same for everyone
        aas_log(aap, false, sprintf('Calculating connectivity between %d voxels.', length(goodVoxelInd)));
                                
        nSess = length(SPM.nscan);
        
        % keep track of correlations and p values
        
        fprintf('Making bix voxel to hold correlations'); % <-- this causes machine to grind to a halt due to using too much RAM
        fconn.rBySession = zeros(length(goodVoxelInd));
       % fconn.pBySession = zeros(length(goodVoxelInd));
        
        pause(10);
        return
       
        % For each session, read in data and do correlations
        for thisSess = 1:nSess
            aas_log(aap, false, sprintf('Processing session %d...', thisSess));
            
            % Figure out which scans belong to this session
            if thisSess==1
                sessIndex = 1:SPM.nscan(thisSess);
            else
                sessIndex = sum(SPM.nscan(1:thisSess-1))+1:sum(SPM.nscan(1:thisSess));
            end
            
            
            % read in resids
            V = spm_vol(strvcat(residNamesImg(sessIndex)));
            [Y, XYZ] = spm_read_vols(V);
            sizeY = size(Y);
                                    
            % Reshape data into time x voxels matrix
            Y = reshape(Y, sizeY(4), prod(sizeY(1:3)));
            
            %...and just use the good voxels
            Y = Y(goodVoxelInd);
            
          
            % [fconn.rBySession(:,:,thisSess), fconn.pBySession(:,:,thisSess)] = corrcoef(Y);
            
            for m = 1:nVoxels
                for n=m:Nvoxels-1
                    [r, p] = corrcoef(Y(:,m), Y(:,n));
                    fconn.rBySession(m,n,thisSess) = r(2);
                end
            end
            
            
        end % looping through sessions
        
        
        % average over sessions for r values
        fconn.r = mean(fconn.rBySession, 3);
        %fconn.r = 
        %fconn. =
        
        
        %% TODO: If ROIwise, read in already-extracted ROIs
        
        
        % output
        [pth, nm, ext] = fileparts(spmName);
        fconnPath = fullfile(pth, 'fconn.mat');
        save(fconnPath, 'fconn');
        
        aap = aas_desc_outputs(aap, subjInd, 'firstlevel_fconn_matrix_mat', fconnPath);
        
        
        
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end




