function [aap,resp]=aamod_roi_extract(aap,task,varargin)

% Extract ROI data from image volume(s), either 4D (e.g., resting state
% timecourses) or 4D (e.g., grey matter signal). Gets data from image files
% themselves; does not require SPM.mat.
%
% INPUT options [defaults]:
%  ROIfile     = path to image containing ROIs
%  ROIvals     = vector of values in ROI image that define ROI [auto]
%  output_raw  = flag to output raw voxel data [0]
%  do_svd      = flag to do SVD calculation [1]
%  NvoxThr     = N (or proportion if <1) non-zero necessary to keep ROI [.5]
%  svd_tol     = tolerance threshold for SVD [1e-6]
%  verbose     = verbose output? [0]
%
% OUTPUT: ROI struct with fields (when requested/apply to data type):
%  ROIval      = value of ROI in mask image
%  XYZ_*       = ROI XYZ coordinates in (*=mask, data) image (mm)
%  XYZcentre_* = centre of ROI XYZ coords in (*=mask, data) image (mm)
%  Nvox_mask   = number of voxels in ROI in mask image
%  Nvox_data   = number of voxels in ROI in data (w/non-zero variance if 4D)
%  rawdata     = raw voxel data (if requested)
%  mean        = mean (timecourse if 4D) over voxels in ROI
%  median      = median (timecourse if 4D) over voxels in ROI
%  svd         = first temporal singular value (if 4D)
%  svd_vox     = first spatial singular value (if 4D)
%  svd_pvar    = % variance accounted for by sv above
%  svd_tol     = svd tolerance for convergence
%
% by Jason Taylor (13 Nov 2012) based on code from Rik Henson
%  + jt (16/Jan/2013): now works for 3D (e.g., GM) data (but still hacky)
%  + rh (05/Feb/2013): now a function
%  + jt (05/Feb/2013): now an aa module
%  + jt (07/Feb/2013): XYZcentre now ROI location in mask (mm)
%  + jt (13/Mar/2013): XYZ and centre now for both mask and data
%  + jt (25/Jul/2013): zero_tol, svd_tol now input options
%  + jt (08/Aug/2013): changed zero_tol to NvoxThr, no_svd to do_svd
%  + jt (22/Oct/2013): per RH, added svd_vox
%  + jt (25/Oct/2013): stop removing voxel means (screws up svd),
%                      multiply SVD by sign of corr with data

resp='';

switch task
    case 'domain'
        resp='session';
        
    case 'description'
        resp='Extract data from ROIs.';
        
    case 'doit'
      
        %% ROIs
        ROIfile    = aap.tasklist.currenttask.settings.ROIfile;
        
        try    ROIvals    = aap.tasklist.currenttask.settings.ROIvals;
        catch, ROIvals = [];
        end
        try    output_raw = aap.tasklist.currenttask.settings.output_raw;
        catch, output_raw = 0;
        end
        try    do_svd = aap.tasklist.currenttask.settings.do_svd;
        catch, do_svd = 1;
        end
        try    NvoxThr = aap.tasklist.currenttask.settings.NvoxThr;
        catch, NvoxThr = 0.5;
        end
        try    NvoxAbsThr = aap.tasklist.currenttask.settings.NvoxAbsThr;
        catch, NvoxAbsThr = 20;
        end
        try    verbose = aap.tasklist.currenttask.settings.verbose;
        catch, verbose = 0;
        end
        try    svd_tol = aap.tasklist.currenttask.settings.svd_tol;
        catch, svd_tol = 1e-6;
        end
        
        % Read ROI volume:
        VV = spm_vol(ROIfile);
        [Y,XYZm] = spm_read_vols(VV);
        ROIfstem = spm_str_manip(ROIfile,'rt');
        
        %% Process
        for in = aap.tasklist.currenttask.inputstreams.stream
            instream = in{1};            
            
            % Get files by stream:
            Datafiles = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,[varargin{:}],instream);
            % Read data:
            VY = spm_vol(Datafiles);
            %VY = [VY{:}];
            
            % Transform ROI XYZ to image space (assumes images coregistered):
            Yinv  = inv(VY(1).mat);
            XYZ   = Yinv(1:3,1:3)*XYZm + repmat(Yinv(1:3,4),1,size(XYZm,2));
            pth = fileparts(Datafiles);
            
            if isempty(ROIvals)
                ROIvals = setdiff(unique(Y(:)),0);
                
                %else
                %if strcmp(ROIvals,'>0');  % could add more special options here...
                %f = find(Y>0);        % but for now, require user to create mask
                %Y(f) = 1;
                %ROIvals = [1];
                %end
                
            end;
            
            Nrois = length(ROIvals);
            
            % SPM progress bar:
            spm_progress_bar('Init',Nrois,'Extracting ROI data','ROIs Complete');
            
            % Gather ROI data:
            ROI = struct();
            for v=1:Nrois,
                
                f = find(Y==ROIvals(v));
                d = spm_get_data(VY,XYZ(:,f)); % Not 100% sure that XYZ is correct (Rik's note)
                Nvox = size(d,2);
                if verbose
                    fprintf('Region %d (%s = %d): %d ',v,ROIfstem,ROIvals(v),length(f));
                end
                ROI(v).ROIval         = ROIvals(v);
                ROI(v).XYZ_mask       = XYZm(:,f);
                ROI(v).XYZcentre_mask = mean(XYZm(:,f),2);
                ROI(v).XYZ_data       = XYZ(:,f);
                ROI(v).XYZcentre_data = mean(XYZ(:,f),2);
                ROI(v).Nvox_mask = Nvox;
                
                % Output raw data?:
                if output_raw
                    ROI(v).rawdata = d;
                end
                
                % Multiple images (e.g., EPI data), extract timecourses:
                % Check for zero-variance voxels:
                
                zero_vox = var(d)==0;
                zero_count = sum(zero_vox);
                ROI(v).Nvox_data = Nvox-zero_count;
                
                % If too many zero-var vox, set ROI data to NaN:
                
                if NvoxThr>0
                    if NvoxThr<1
                        nvox_crit = ceil(NvoxThr*Nvox); % proportion
                    else
                        nvox_crit = NvoxThr;            % absolute N
                    end
                else
                    nvox_crit = 0;
                end
                ROI(v).NvoxThr   = NvoxThr;
                ROI(v).Nvox_crit = nvox_crit;
                
                if (ROI(v).Nvox_data < nvox_crit) || (ROI(v).Nvox_data < NvoxAbsThr)
                    if verbose
                        fprintf('(%d nonzero) voxels -- FAILED (<%d)!\n',ROI(v).Nvox_data,nvox_crit);
                    end
                    ROI(v).mean     = NaN(size(d,1),1);
                    ROI(v).median   = NaN(size(d,1),1);
                    ROI(v).svd      = NaN(size(d,1),1);
                    ROI(v).svd_vox  = NaN(size(d,2),1);
                    ROI(v).svd_pvar = NaN;
                    ROI(v).svd_tol  = svd_tol;
                else
                    
                    % Remove zero-variance voxels:
                    f = setdiff(f,zero_vox);
                    d = spm_get_data(VY,XYZ(:,f));
                    
                    if verbose
                        fprintf('(%d nonzero) voxels\n',Nvox-zero_count);
                    end
                    
                    % Remove voxel-wise mean:
                    %dd = d-repmat(mean(d,1),size(d,1),1);
                    % commented out, replaced dd with d below jt 25/Oct/2013
                    
                    % MEAN/MEDIAN:
                    ROI(v).mean   = mean(d,2);
                    ROI(v).median = median(d,2);
                    
                    % SVD:
                    if do_svd
                        [U,S,V] = spm_svd(d,svd_tol);
                        
                        if isempty(S)
                            if verbose
                                fprintf('..SVD FAILED!\n');
                            end
                            ROI(v).svd      = NaN(size(d,1),1);
                            ROI(v).svd_vox  = NaN(size(d,2),1); % per RH
                            ROI(v).svd_pvar = NaN;
                            ROI(v).svd_tol  = svd_tol;
                        else
                            % Get 1st temporal, spatial modes:
                            U = full(U(:,1));
                            V = full(V(:,1));
                            % Change sign to sign of corr w/data:
                            cor = corr(U,ROI(v).mean);
                            U = sign(cor)*U;
                            V = sign(cor)*V;
                            % Store 'em:
                            ROI(v).svd      = U;
                            ROI(v).svd_vox  = V;
                            ROI(v).svd_pvar = S(1)/sum(full(S(:)));
                            ROI(v).svd_tol  = svd_tol;
                        end
                        %if isempty(S)
                        %    svd_tol = 1e-9; % could increase tolerance?
                        %    [U,S,V] = spm_svd(d,svd_tol);
                        %end
                    else
                        ROI(v).svd      = NaN(size(d,1),1);
                        ROI(v).svd_vox  = NaN(size(d,2),1); % per RH
                        ROI(v).svd_pvar = NaN;
                        ROI(v).svd_tol  = NaN;
                    end
                end
                
                % Update SPM progress bar (every 5%):
                if mod(round(100*v/Nrois),5)==0
                    spm_progress_bar('Set',v);
                end
                
                
            end
            
            % Clear SPM progress bar:
            spm_progress_bar('Clear');
            
            % Describe outputs:
            outstream = ['roidata_' instream];
            outfile = fullfile(pth,['ROI_' instream '.mat']);
            save(outfile,'ROI');
            aap = aas_desc_outputs(aap,aap.tasklist.currenttask.domain,[varargin{:}],outstream,outfile);
            
        end
end


end

