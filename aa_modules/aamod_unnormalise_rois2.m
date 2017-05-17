% Get diffusion data streams (dti_FA).
% unnormalise the seed or target (MNI) specified in the .xml file
% to the individual diffusion space (dti_FA)

function [aap resp]=aamod_unnormalise_rois2(aap,task,subjind,diffsessind)

resp='';

switch task
    case 'report'
    case 'doit'
        fslext=aas_getfslext(aap);
        dsesspth= aas_getpath_bydomain(aap,'diffusion_session',[subjind diffsessind]);
        roipth=aap.tasklist.currenttask.settings.roi_path;
        
        [aap FAref]=aas_gunzip(aap,aas_getfiles_bystream(aap,'diffusion_session',[subjind diffsessind],'dti_FA'));
        invNormPar=aas_getfiles_bystream(aap,'subject',subjind,'normalisation_seg_inv_sn');
        
        rois=aap.tasklist.currenttask.settings.rois;
        
        % Set up parameters for inverse normalisation
        %   Infinite bounding box needed to avoid ROIs being cut...
        aap.spm.defaults.normalise.write.bb = [Inf Inf Inf; Inf Inf Inf];
        aap.spm.defaults.normalise.write.vox = [Inf Inf Inf];
        
        %   Get realignment defaults
        defs = aap.spm.defaults.realign;
        
        %   Flags to pass to routine to create resliced images
        %   (spm_reslice)
        resFlags = struct(...
            'interp', defs.write.interp,...       % interpolation type
            'wrap', defs.write.wrap,...           % wrapping info (ignore...)
            'mask', defs.write.mask,...           % masking (see spm_reslice)
            'which', 1,...     % what images to reslice
            'mean', 0);           % write mean image
        
        
        % For each  ROI...
        roiind=1;
        for roigroupind=1:length(rois.roi)
            for roirangeind=rois.roi(roigroupind).low:rois.roi(roigroupind).high
                roifn=sprintf(rois.roi(roigroupind).filename,roirangeind);
                [pth nme ext]=aas_fileparts(roifn);
                
                % First transform
                outfn{roiind}=fullfile(dsesspth, ['w' nme '.nii']);
                outfn_r{roiind}=fullfile(dsesspth, ['rw' nme '.nii']);
                
                % Name of file
                outfn_mni{roiind}=fullfile(dsesspth, [nme ext]);
                
                % Copy across roi
                copyfile(fullfile(roipth,roifn),outfn_mni{roiind})
                
                % Unpack if necessary
                [aap fn]=aas_gunzip(aap,outfn_mni{roiind});
                outfn_mni{roiind}=fn;
                
                
                % Inverse normaliase
                % Un-Normalise ROI to structural initially...
                spm_write_sn(outfn_mni{roiind}, invNormPar, aap.spm.defaults.normalise.write);
                
                % Reslice to FA...
                spm_reslice(strvcat(FAref, outfn{roiind}), resFlags)
                
                % And binarize
                outfn_bin{roiind}=fullfile(dsesspth, [nme '_diffusion_space_bin' fslext]);
                cmd=sprintf('fslmaths %s -thr 0.50 -bin %s', outfn_r{roiind},outfn_bin{roiind});
                [s w]=aas_runfslcommand(aap,cmd);
                if (s)
                    aas_log(aap,true,sprintf('Error executing\n  %s\nof\n%s',cmd,w));
                end;
                roiind=roiind+1;
            end;
            
        end;
        %
        aap=aas_desc_outputs(aap,'diffusion_session',[subjind diffsessind], aap.tasklist.currenttask.outputstreams.stream{1}, outfn_bin);
        aap=aas_desc_outputs(aap,'study',[], [aap.tasklist.currenttask.outputstreams.stream{1} '_mni'], outfn_mni);
        
end
end

