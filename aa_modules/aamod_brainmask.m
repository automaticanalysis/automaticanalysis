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
        
        % In cases where the input stream isn't a subject domain (e.g.,
        % meanEPI), try looking in the first session. hacky!
		
        try
            gmImg = aas_getfiles_bystream(aap, subjInd, inStream);
        catch
            gmImg = aas_getfiles_bystream(aap, subjInd, 1, inStream);
        end
        
        % Get path for structural image, and make path for thresholded
		
        [structPath, structName, structExt] = fileparts(gmImg);
        threshFilename = fullfile(structPath, sprintf('binary_%s%s', structName, structExt));
        
        V = spm_vol(gmImg);
        [Y, XYZ] = spm_read_vols(V);
        
        % Binarize based on threshold and rename
		
        if ~isempty(thresh)
            Y = Y > thresh;
        else
            Y = Y > mean(Y(:));
        end
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
				
				% assume 'reslice' = a stream name
				%
				% In cases where the stream isn't a subject domain (e.g.,
				% meanEPI), try looking in the first session. hacky!

				try
					resliceImg = aas_getfiles_bystream(aap, subjInd, reslice);
				catch
					aas_log(aap,false,'**** retrying file read in session 1... ');
					resliceImg = aas_getfiles_bystream(aap, subjInd, 1, reslice);
				end
				
				% if ref image is multi-volume, use the first file
			   
				if (iscell(resliceImg)) resliceImg = resliceImg{1}; end
				
				% if the ref image is 4D, use the first vol
				% (this should be harmless if not 4D)
				
				resliceImg =  [resliceImg ',1'];
				
			end
           
			P = strvcat(resliceImg, threshFilename);    
			spm_reslice(P, resliceOpts);
			spm_unlink(threshFilename);   
			threshFilename = fullfile(structPath, sprintf('rbinary_%s%s', structName, structExt));
			
		end
        
        % If requested, smooth outname and re-threshold > 0
		
        if ~isempty(fwhm) && fwhm>0
            [threshPath, threshName, threshExt] = fileparts(threshFilename);
            smoothImg = fullfile(threshPath, ['s' threshName threshExt]);
            spm_smooth(threshFilename, smoothImg, fwhm);
            % re-threshold
            Vsmooth = spm_vol(smoothImg);
            [Ysmooth,~] = spm_read_vols(Vsmooth);
            Ysmooth = Ysmooth > 0;
            Vsmooth.fname = threshFilename;
            spm_write_vol(Vsmooth, Ysmooth);
			spm_unlink(smoothImg);
        end

        aap = aas_desc_outputs(aap, subjInd, outStream, threshFilename);
        
    case 'checkrequirements'
        
end
