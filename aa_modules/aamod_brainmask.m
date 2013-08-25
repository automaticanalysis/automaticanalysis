% AA module - Create a binary brainmask based on segmentation parameters
% [aap,resp]=aamod_brainmask(aap, task, subj)
%
% 
function [aap,resp]=aamod_brainmask(aap, task, subjInd)
resp='';

switch task
    case 'report'
        
    case 'doit'
        thresh = aap.tasklist.currenttask.settings.thresh; % how to threshold image
        reslice = aap.tasklist.currenttask.settings.reslice; % if not empty, reslice to match this image using NN (to keep binary)
        fwhm = aap.tasklist.currenttask.settings.fwhm; % if not empty, smooth and re-threshold > 0 (cheat way of dilating mask)
        
        % Input stream could be e.g. normalized to MNI or subject-specific,
        % which have different stream names
        inStream = aap.tasklist.currenttask.inputstreams(1).stream{1};
        outStream = aap.tasklist.currenttask.outputstreams(1).stream{1};
        
        gmImg = aas_getfiles_bystream(aap, subjInd, inStream);
        
        % Get path for structural image, and make path for thresholded
        [structPath, structName, structExt] = fileparts(gmImg);
        threshFilename = fullfile(structPath, sprintf('binary_%s%s', structName, structExt));
        
        
        V = spm_vol(gmImg);
        [Y, XYZ] = spm_read_vols(V);
        
        % Binarize based on threshold and rename
        Y = Y > thresh;
        V.fname = threshFilename;
        V.descrip = sprintf('%s thresholded at %g by aamod_brainmask', gmImg, thresh);
        
        % Write out un-resliced image
        spm_write_vol(V,Y);
        
        
        % If reslicing specified, do that. Depending on the XML file, this
        % could be to a named file, or to an AA stream. Trying to write
        % code below to handle either somewhat gracefully.
        if ~isempty(reslice)
            
            resliceOpts = [];
            resliceOpts.mask = false;
            resliceOpts.mean = false;
            resliceOpts.interp = 0; % NN
            resliceOpts.which = 1; % don't reslice the first image
            resliceOpts.wrap = [0 0 0];
            resliceOpts.prefix = 'r';
            
            
           if exist(reslice, 'file')
               resliceImg = reslice;               
           else
               % assume stream
               resliceImg = aas_getfiles_bystream(aap, subjInd, reslice); 
               resliceImg = resliceImg{1}; % use the first one if more than one (e.g. epis)
           end
           
           P = strvcat(reslice, threshFilename);
               
           spm_reslice(P, resliceOpts);
           
           threshFilename = fullfile(structPath, sprintf('rbinary_%s%s', structName, structExt));
        end
        
        
        % If requested, smooth outname and re-threshold > 0
        if ~isempty(fwhm) && fwhm>0
            [threshPath, threshName, threshExt] = fileparts(threshFilename);
            smoothImg = fullfile(threshPath, ['s' threshName threshExt]);
            spm_smooth(threshFilename, smoothImg, fwhm);
            
            % re-threshold
            Vsmooth = spm_vol(smoothImg);
            [Ysmooth, XYZsmooth] = spm_read_vols(Vsmooth);
            Ysmooth = Ysmooth > 0;
            Vsmooth.fname = threshFilename;
            spm_write_vol(Vsmooth, Ysmooth);            
        end
        
        
        aap = aas_desc_outputs(aap, subjInd, outStream, threshFilename);
        
    case 'checkrequirements'
        
end
