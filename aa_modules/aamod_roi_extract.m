function [aap,resp]=aamod_roi_extract(aap,task,varargin)

% Extract ROI data from image volume(s), either 4D (e.g., resting state
% timecourses) or 4D (e.g., grey matter signal). Gets data from image files
% themselves; does not require SPM.mat.
%
% OUTPUT: ROI struct with fields (when requested/apply to data type):
%  ROIfile     = ROI file
%  ROIval      = value of ROI in mask image
%  XYZ         = ROI XYZ coordinates in mm (mask or data depending on "mask_space")
%  XYZcentre   = centre of ROI XYZ coords in mm (mask or data depending on "mask_space")
%  Nvox        = number of voxels in ROI (in mask or data depending on "mask_space")
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
%  + jc (27/Jul/2018): refactor to get ROIfile from first input stream instead
%                      of task settings

resp='';

switch task
    case 'domain'
        resp='session';
        
    case 'description'
        resp='Extract data from ROIs.';
        
    case 'doit'
      
        %% Settings
        ROIvals    = aap.tasklist.currenttask.settings.ROIvals;
        
        mask_space = aap.tasklist.currenttask.settings.mask_space;
        output_raw = aap.tasklist.currenttask.settings.output_raw;        
        
        do_svd = aap.tasklist.currenttask.settings.do_svd;
        NvoxThr = aap.tasklist.currenttask.settings.NvoxThr;
        NvoxAbsThr = aap.tasklist.currenttask.settings.NvoxAbsThr;
        svd_tol = aap.tasklist.currenttask.settings.svd_tol;
        
        %% Process
        inputstreams = aas_getstreams(aap,'input');
        % first input stream is the roi stream, all others are data
        [roistream, inputstreams] = deal(inputstreams(1),inputstreams(2:end));
        for in = inputstreams
            instream = in{1};
            
            % Get data (assumes all data files are aligned)
            VY = spm_vol(aas_getfiles_bystream_multilevel(aap,aap.tasklist.currenttask.domain,[varargin{:}],instream));
            Yinv  = inv(VY(1).mat);
            [~, yXYZmm] = spm_read_vols(VY(1));
            
            % Get ROI
            ROI = [];
            for fn = cellstr(aas_getfiles_bystream_multilevel(aap,aap.tasklist.currenttask.domain,[varargin{:}],roistream{1}))'

                ROIfstem = spm_str_manip(fn{1},'rt');
                VM = spm_vol(fn{1});
                Minv = inv(VM.mat);
                [YM,mXYZmm] = spm_read_vols(VM);
                
                % Transform ROI XYZ in mm to voxel indices in data:
                yXYZind = Yinv(1:3,1:3)*mXYZmm + repmat(Yinv(1:3,4),1,size(mXYZmm,2));
                % Transform data XYZ in mm to voxel indices in mask:
                mXYZind = Minv(1:3,1:3)*yXYZmm + repmat(Minv(1:3,4),1,size(yXYZmm,2));
                % Transform data XYZ in mm to voxel indices in data:
                yyXYZind = Yinv(1:3,1:3)*yXYZmm + repmat(Yinv(1:3,4),1,size(yXYZmm,2));
                
                yYM = spm_get_data(VM,mXYZind);
                
                if mask_space
                    % YM = YM;
                    XYZmm = mXYZmm;
                    XYZind = yXYZind;
                else
                    YM = yYM;
                    XYZmm = yXYZmm;
                    XYZind = yyXYZind;
                end
                
                if isempty(ROIvals)
                    ROIvals = setdiff(unique(YM(:)),0);
                    ROIvals(isnan(ROIvals)) = [];
                    %else
                    %if strcmp(ROIvals,'>0');  % could add more special options here...
                    %f = find(Y>0);        % but for now, require user to create mask
                    %Y(f) = 1;
                    %ROIvals = [1];
                    %end
                end
                Nrois = length(ROIvals);
                
                % SPM progress bar:
                spm_progress_bar('Init',Nrois,'Extracting ROI data','ROIs Complete');
                
                % Gather ROI data:
                n0 = numel(ROI);
                for nr = 1:Nrois
                    nrAll = n0+nr; % append
                    ROI(nrAll).ROIfile  = ROIfstem;
                    ROI(nrAll).ROIval   = ROIvals(nr);
                    
                    f = find(YM==ROIvals(nr));
                    d = spm_get_data(VY,XYZind(:,f)); % Think this is correct!
                    ROI(nrAll).XYZ       = XYZmm(:,f);
                    ROI(nrAll).XYZcentre = mean(XYZmm(:,f),2);
                    
                    aas_log(aap,false,sprintf('INFO: Region %d (%s = %d): %d ',nr,ROIfstem,ROIvals(nr),length(f)));
                    
                    % Output raw data?:
                    if output_raw
                        ROI(nrAll).rawdata = d;
                    end
                    
                    % Multiple images (e.g., EPI data), extract timecourses:
                    % Check for zero-variance voxels:
                    Nvox = size(d,2);
                    zero_vox = (var(d)==0 | isnan(var(d)));
                    zero_count = sum(zero_vox);
                    ROI(nrAll).Nvox = Nvox-zero_count;
                    
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
                    ROI(nrAll).NvoxThr   = NvoxThr;
                    ROI(nrAll).Nvox_crit = nvox_crit;
                    
                    % Defaults
                    ROI(nrAll).mean     = NaN(size(d,1),1);
                    ROI(nrAll).median   = NaN(size(d,1),1);
                    ROI(nrAll).svd      = NaN(size(d,1),1);
                    ROI(nrAll).svd_vox  = NaN(size(d,2),1);
                    ROI(nrAll).svd_pvar = NaN;
                    ROI(nrAll).svd_tol  = NaN;
                    
                    if (ROI(nrAll).Nvox < nvox_crit) || (ROI(nrAll).Nvox < NvoxAbsThr)
                        aas_log(aap,false,sprintf('INFO: (%d nonzero) voxels -- FAILED (<%d or %d)!\n',ROI(nrAll).Nvox,nvox_crit,NvoxAbsThr));
                    else
                        
                        % Remove zero-variance voxels:
                        f(zero_vox) = [];
                        d = spm_get_data(VY,XYZind(:,f));
                        
                        aas_log(aap,false,sprintf('INFO: (%d nonzero) voxels\n',Nvox-zero_count));
                        
                        % MEAN/MEDIAN:
                        ROI(nrAll).mean   = mean(d,2);
                        ROI(nrAll).median = median(d,2);
                        
                        % SVD:
                        if do_svd
                            ROI(nrAll).svd_tol  = svd_tol;
                            [U,S,V] = spm_svd(d,svd_tol);
                            
                            if isempty(S)
                                aas_log(aap,false,sprintf('WARNING: ..SVD FAILED!\n'));
                            else
                                % Get 1st temporal, spatial modes:
                                U = full(U(:,1));
                                V = full(V(:,1));
                                % Change sign to sign of corr w/data:
                                cor = corr(U,ROI(nrAll).mean);
                                U = sign(cor)*U;
                                V = sign(cor)*V;
                                % Store 'em:
                                ROI(nrAll).svd      = U;
                                ROI(nrAll).svd_vox  = V;
                                ROI(nrAll).svd_pvar = S(1)/sum(full(S(:)));
                                ROI(nrAll).svd_tol  = svd_tol;
                            end
                            %if isempty(S)
                            %    svd_tol = 1e-9; % could increase tolerance?
                            %    [U,S,V] = spm_svd(d,svd_tol);
                            %end
                        end
                    end
                    
                    % Update SPM progress bar (every 5%):
                    if mod(round(100*nr/Nrois),5)==0
                        spm_progress_bar('Set',nr);
                    end
                    
                    
                end
                
                % Clear SPM progress bar:
                spm_progress_bar('Clear');
            end
            % Describe outputs:
            instream = strsplit(instream,'.'); instream = instream{end};
            outstream = ['roidata_' instream];
            outfile = fullfile(aas_getpath_bydomain(aap,aap.tasklist.currenttask.domain,[varargin{:}]),['ROI_' instream '.mat']);
            save(outfile,'ROI');
            aap = aas_desc_outputs(aap,aap.tasklist.currenttask.domain,[varargin{:}],outstream,outfile);

        end

    case 'checkrequirements'
        inputstreams = aas_getstreams(aap,'input');
        % first input stream is the roi stream, all others are data
        inputstreams = inputstreams(2:end);
        out = aas_getstreams(aap,'output'); if ~iscell(out), out = {out}; end
        for s = 1:numel(inputstreams)
            instream = strsplit(inputstreams{s},'.'); instream = instream{end};
            if ~strcmp(out{s},['roidata_' instream])
                aap = aas_renamestream(aap,...
                    aap.tasklist.currenttask.name,out{s},...
                    ['roidata_' instream],'output');
                aas_log(aap,false,['INFO: ' aap.tasklist.currenttask.name ...
                    ' output stream: ''roidata_' instream '''']);
            end            
        end        
end


end

