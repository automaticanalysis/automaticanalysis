% AA module - UNNORMALISE_ROI Unnormalise a particular ROI for a given subject
%
% Alejandro Vicente-Grabovetsky Dec-2008

function [aap,resp] = aamod_reslice_rois(aap, task, subj)

resp='';

switch task
    case 'doit'
        
        % Get meanEPImage
        mEPIimg = aas_getfiles_bystream(aap,subj,'meanepi');
        if size(mEPIimg,1) > 1
            aas_log(aap, false, 'Found more than 1 mean functional images, using first.');
        end
        mEPIimg = deblank(mEPIimg(1,:));
        
        % Get roi image(s)
        Rimg = aas_getfiles_bystream(aap,subj,'roi');
        
        % Test format, ensure it is V.dt = [2 0];
        for c=1:size(Rimg,1)
            img2mask(Rimg(c,:));
        end
        
        % Get realignment defaults
        defs = aap.spm.defaults.realign;
        
        % Flags to pass to routine to create resliced images
        % (spm_reslice)
        resFlags = struct(...
            'interp', defs.write.interp,...       % interpolation type
            'wrap', defs.write.wrap,...           % wrapping info (ignore...)
            'mask', 0, ... % defs.write.mask,...           % masking (see spm_reslice)
            'which', 1,...     % what images to reslice
            'mean', 0);           % write mean image
        
        outstream = {};
        
        % Loop through all rois
        for r = 1:length(Rimg(:,1))
            Rcurr = deblank(Rimg(r,:));
            [Rpth, Rfn, Rext] = fileparts(Rcurr);
                                                 
            % Reslice to EPI...
            spm_reslice({mEPIimg, Rcurr}, resFlags);            
            Rcurr = fullfile(Rpth, ['r' Rfn, Rext]);
            
            % Check the content of the resliced ROI!
            V = spm_vol(Rcurr);
            Y = spm_read_vols(V);
            if nansum(Y(:)>0) == 0
                aas_log(aap,false,'WARNING: This ROI contains no voxels')
            end
            
            outstream = [outstream, Rcurr];
        end
        
        aap=aas_desc_outputs(aap,subj,'roi',outstream);
        
        %% Save graphical output to common diagnostics directory
        subjname = aas_prepare_diagnostic(aap,subj);
        OVERcolours = aas_colours;
        
        spm_check_registration(mEPIimg)
        
        % Add segmentations...
        for t = 1:length(outstream)
            spm_orthviews('addcolouredimage',1,outstream{t}, OVERcolours{t})
        end
        
        spm_orthviews('reposition', [0 0 0])
        print('-djpeg','-r150',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' subjname '.jpeg']));        
end