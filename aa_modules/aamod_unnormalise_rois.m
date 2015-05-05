% AA module - UNNORMALISE_ROI Unnormalise a particular ROI for a given subject
%
% Alejandro Vicente-Grabovetsky Dec-2008

function [aap,resp] = aamod_unnormalise_rois( aap, task, subj)

resp='';

switch task
    case 'checkrequirements'
        ROIlist{1} = aap.tasklist.currenttask.settings.ROIlist;
        if ~exist(ROIlist{1},'file')
            ROIs = textscan(ROIlist{1},'%s','delimiter',','); ROIs = ROIs{1};
            ROIlist = cell(numel(ROIs),1);
            for r = 1:numel(ROIs)
                ROIlist{r} = fullfile(aap.directory_conventions.ROIdir, ROIs{r});
            end
        end
        
        for r = 1:numel(ROIlist)
            if ~exist(ROIlist{r}, 'file')
                aas_log(aap,true,sprintf('The ROI %s file does not exist', ...
                    ROIlist{r}))
            end
            [junk, fn, ext] = fileparts(ROIlist{r});
            if strcmp(ext,'.img')
                aas_log(aap,true,sprintf('ROIs of analyze format are not supported - convert following to .nii\n  %s',ROIlist{r}));
            end;
        end;
    case 'doit'
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
                    aas_log(aap,false,'No segmentation mask will be applied...\n');
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
        
        mEPIimg = aas_getfiles_bystream(aap,subj,'meanepi');
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
        normFlags.interp = aap.tasklist.currenttask.settings.interp;
        normFlags.bb = [Inf Inf Inf; Inf Inf Inf];
        normFlags.vox = [Inf Inf Inf];
        
        % Reslice flags
        resFlags = aap.spm.defaults.coreg.write;
        resFlags.interp = aap.tasklist.currenttask.settings.interp;
        resFlags.mask = 0;
        resFlags.which = [1 0];
        
        % Get the ROIs from .xml
        ROIlist{1} = aap.tasklist.currenttask.settings.ROIlist;
        if ~exist(ROIlist{1},'file')
            ROIs = textscan(ROIlist{1},'%s','delimiter',','); ROIs = ROIs{1};
            ROIlist = cell(numel(ROIs),1);
            for r = 1:numel(ROIs)
                ROIlist{r} = fullfile(aap.directory_conventions.ROIdir, ROIs{r});
            end
        end
        outstream = '';
        
        % Loop through all ROIs
        for r = 1:numel(ROIlist)
            if ~exist(ROIlist{r}, 'file')
                aas_log(aap,true,sprintf('The ROI %s file does not exist', ...
                    ROIlist{r}))
            end
            % Copy to structural dir...
            unix(['cp ' ROIlist{r} ' ' fullfile(aas_getsubjpath(aap,subj), aap.directory_conventions.structdirname) ]);
            
            [junk, fn, ext] = fileparts(ROIlist{r});
            % New location for ROI...
            roi_fn = fullfile(aas_getsubjpath(aap,subj), aap.directory_conventions.structdirname, [fn, ext]);
            
            % Un-Normalise ROI to structural initially...
            spm_write_sn(roi_fn, invNormPar, normFlags);
            % Delete ROI in MNI space...
            unix(['rm -rf ' roi_fn]);
            roi_fn = fullfile(aas_getsubjpath(aap,subj), aap.directory_conventions.structdirname, ['w' fn, ext]);
            
            % Reslice to EPI...
            spm_reslice(strvcat(mEPIimg, roi_fn), resFlags)
            % Delete ROI in MNI space...
            unix(['rm -rf ' roi_fn]);
            roi_fn = fullfile(aas_getsubjpath(aap,subj), aap.directory_conventions.structdirname, ['rw' fn, ext]);
            
            if ~isempty(SEGimg)
                % Now mask our ROIs by segmented mask (t for trimmed)
                conjMask(roi_fn, Mimg, [0.01 0.01], 't');
                unix(['rm -rf ' roi_fn]);
                roi_fn = fullfile(aas_getsubjpath(aap,subj), aap.directory_conventions.structdirname, ['trw' fn, ext]);
            end
            
            V = spm_vol(roi_fn);
            Y = spm_read_vols(V);
            if nansum(Y(:)>0) == 0
                aas_log(aap,false,'This ROI contains no voxels')
            end
            
            outstream = strvcat(outstream, roi_fn);
        end
        
        aap=aas_desc_outputs(aap,subj,'rois',outstream);
        % Diag
        if strcmp(aap.options.wheretoprocess,'localsingle')
            aas_checkreg(aap,subj,'rois','structural');
        end
    case 'report'
        localpath = aas_getpath_bydomain(aap,'subject',subj);
        d = dir(fullfile(localpath,'diagnostic_aas_checkreg_*'));
        if isempty(d)
            struct = aas_getfiles_bystream_dep(aap,'subject',subj,'structural');
            aas_checkreg(aap,subj,'rois',struct);
        end
        fdiag = dir(fullfile(localpath,'diagnostic_aas_checkreg_*.jpg'));
        for d = 1:numel(fdiag)
            aap = aas_report_add(aap,subj,'<table><tr><td>');
            imgpath = fullfile(localpath,fdiag(d).name);
            aap=aas_report_addimage(aap,subj,imgpath);
            [p, f] = fileparts(imgpath); avipath = fullfile(p,[strrep(f(1:end-2),'slices','avi') '.avi']);
            if exist(avipath,'file'), aap=aas_report_addimage(aap,subj,avipath); end
            aap = aas_report_add(aap,subj,'</td></tr></table>');
        end
end