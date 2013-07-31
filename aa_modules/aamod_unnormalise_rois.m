% AA module - UNNORMALISE_ROI Unnormalise a particular ROI for a given subject
%
% Alejandro Vicente-Grabovetsky Dec-2008

function [aap,resp] = aamod_unnormalise_rois( aap, task, subj)

resp='';

switch task
    case 'checkrequirements'
        ROIlist=getROIlist(aap);
        for r = 1:length(ROIlist)
            if ~exist(ROIlist{r}, 'file')
                error(sprintf('The ROI %s file does not exist', ...
                    ROIlist{r}))
            end
            [junk, fn, ext] = fileparts(ROIlist{r});
            if strcmp(ext,'.img')
                aas_log(aap,true,sprintf('ROIs of analyze format are not supported - convert following to .nii\n  %s',ROIlist{r}));
            end;
        end;
    case 'doit'
        
        % Structural directory...
        if ~exist(fullfile(aas_getsubjpath(aap,subj), 'structurals'), 'dir')
            mkdir(fullfile(aas_getsubjpath(aap,subj), 'structurals'))
        end
        
        % Get segmentation masks we wish to use
        try
            SEGimg = aas_getfiles_bystream(aap,subj,'segmasksExclusive');
        catch
            try
                SEGimg = aas_getfiles_bystream(aap,subj,'segmasksStrict');
            catch
                try
                    SEGimg = aas_getfiles_bystream(aap,subj,'segmasksZero');
                catch
                    SEGimg = [];
                    fprintf('No segmentation mask will be applied...\n');
                end
            end
        end
        
        if ~isempty(SEGimg)
            % Select the Segmentation mask we wish to reslice to/use...
            for m = 1:size(SEGimg,1)
                if strfind(SEGimg(m,:), sprintf('rc%d', aap.tasklist.currenttask.settings.maskNum))
                    Mimg = deblank(SEGimg(m,:));
                    break
                end
            end
        end
        
        mEPIimg = aas_getfiles_bystream(aap,subj,1,'meanepi');
        if isempty(mEPIimg)
            aas_log(aap, true, 'Problem finding mean functional image.');
        elseif size(mEPIimg,1) > 1
            aas_log(aap, false, 'Found more than 1 mean functional images, using first.');
        end
        mEPIimg = deblank(mEPIimg(1,:));
        
        % Get inverse normalisation parameters
        invNormPar = aas_getfiles_bystream(aap,subj,'normalisation_seg_inv_sn');
        % Normalisation flags...
        % Infinite bounding box needed to avoid ROIs being cut...
        aap.spm.defaults.normalise.write.bb = [Inf Inf Inf; Inf Inf Inf];
        aap.spm.defaults.normalise.write.vox = [Inf Inf Inf];
        
        % Get realignment defaults
        defs = aap.spm.defaults.realign;
        
        % Flags to pass to routine to create resliced images
        % (spm_reslice)
        resFlags = struct(...
            'interp', defs.write.interp,...       % interpolation type
            'wrap', defs.write.wrap,...           % wrapping info (ignore...)
            'mask', defs.write.mask,...           % masking (see spm_reslice)
            'which', 1,...     % what images to reslice
            'mean', 0);           % write mean image
        
        % [AVG] Let's try without extraparameters, so we can have things in the
        % tasklist instead...
        % EP = aap.tasklist.currenttask.extraparameters;
        
        % Get the ROIs from .xml
        ROIstr = aap.tasklist.currenttask.settings.ROIlist;
        ROIlist = {};
        while ~isempty(ROIstr)
            [tmp, ROIstr] = strtok(ROIstr,',');
            ROIlist = [ROIlist fullfile(...
                aap.directory_conventions.ROIdir, tmp)];
        end
        
        outstream = '';
        
        % Loop through all ROIs
        for r = 1:length(ROIlist)
            if ~exist(ROIlist{r}, 'file')
                error(sprintf('The ROI %s file does not exist', ...
                    ROIlist{r}))
            end
            % Copy to structural dir...
            unix(['cp ' ROIlist{r} ' ' fullfile(aas_getsubjpath(aap,subj), 'structurals') ]);
            
            [junk, fn, ext] = fileparts(ROIlist{r});
            % New location for ROI...
            roi_fn = fullfile(aas_getsubjpath(aap,subj), 'structurals', [fn, ext]);
            
            % Un-Normalise ROI to structural initially...
            spm_write_sn(roi_fn, invNormPar, aap.spm.defaults.normalise.write);
            % Delete ROI in MNI space...
            unix(['rm -rf ' roi_fn]);
            roi_fn = fullfile(aas_getsubjpath(aap,subj), 'structurals', ['w' fn, ext]);
            
            % Reslice to EPI...
            spm_reslice(strvcat(mEPIimg, roi_fn), resFlags)
            % Delete ROI in MNI space...
            unix(['rm -rf ' roi_fn]);
            roi_fn = fullfile(aas_getsubjpath(aap,subj), 'structurals', ['rw' fn, ext]);
            
            if ~isempty(SEGimg)
                % Now mask our ROIs by segmented mask (t for trimmed)
                conjMask(roi_fn, Mimg, [0 0.01], 't');
                unix(['rm -rf ' roi_fn]);
                roi_fn = fullfile(aas_getsubjpath(aap,subj), 'structurals', ['trw' fn, ext]);
            end
            
            V = spm_vol(roi_fn);
            Y = spm_read_vols(V);
            if nansum(Y(:)>0) == 0
                warning('This ROI contains no voxels')
            end
            
            outstream = strvcat(outstream, roi_fn);
        end
        
        % Diagnostic image?
        % Save graphical output to common diagnostics directory
        if ~exist(fullfile(aap.acq_details.root, 'diagnostics'), 'dir')
            mkdir(fullfile(aap.acq_details.root, 'diagnostics'))
        end
        mriname = strtok(aap.acq_details.subjects(subj).mriname, '/');
        try
            %% Draw mean EPI...
            spm_check_registration(mEPIimg)
            
            % This will only work for 1-7 ROIs
            OVERcolours = {[1 0 0], [0 1 0], [0 0 1], ...
                [1 1 0], [1 0 1], [0 1 1], [1 1 1]};
            
            % Add un-normalised ROIs...
            for r = 1:size(outstream,1)
                spm_orthviews('addcolouredimage',1,outstream(r,:), OVERcolours{r})
            end
            %% Diagnostic VIDEO of segmentations
            aas_checkreg_avi(aap, subj, 2)
            
            spm_orthviews('reposition', [0 0 0])
            
            try figure(spm_figure('FindWin', 'Graphics')); catch; figure(1); end;
            set(gcf,'PaperPositionMode','auto')
            print('-djpeg','-r75',fullfile(aap.acq_details.root, 'diagnostics', ...
                [mfilename '__' mriname '.jpeg']));
        catch
        end
        
        aap=aas_desc_outputs(aap,subj,'rois',outstream);
end