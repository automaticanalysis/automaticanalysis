% aa module - Mask images from a segmented structural
% 1) Take segmentation images
% 2) Threshold them at one of the 3 different levels (zero level, strict [e.g. 99%],
% or an exclusive [highest tissue probability wins])

function [aap,resp] = aamod_mask_fromsegment(aap,task,subj)

resp = '';

switch task
    case 'doit'
        
        inStreams = aas_getstreams(aap,'input');
        refimg = aas_getfiles_bystream(aap,subj,inStreams{1});
        segimg = char(...
            aas_getfiles_bystream(aap,subj,inStreams{2}),...
            aas_getfiles_bystream(aap,subj,inStreams{3}),...
            aas_getfiles_bystream(aap,subj,inStreams{4})...
            );
        
        
        %% 1) RESLICE
        % Get defaults
        resFlags = aap.spm.defaults.coreg.write;
        resFlags.which = 1;
        resFlags.mean = 0;
        resFlags.mask = 1;
        
        spm_reslice(strvcat(refimg,segimg),resFlags);
        
        %% 2) THRESHOLD
        
        %% IN PROGRESS
        outstream = '';
        
        % Load the correct resliced file!
        V = cell(1, size(segimg,1));
        Y = cell(1, size(segimg,1));
        for a = 1:size(segimg,1)
            V{a} = spm_vol(spm_file(segimg(a,:),'prefix',resFlags.prefix));
            Y{a} = spm_read_vols(V{a});
        end
        
        switch aap.tasklist.currenttask.settings.threshold
            case 'zero'
                %% A) Zero thresholding is easy...
                for a = 1:size(segimg,1)
                    Y{a} = Y{a} > 0;

                    V{a}.fname = spm_file(segimg(a,:),'prefix','Z_r');
                    outstream = strvcat(outstream, V{a}.fname); % Save to stream...
                    spm_write_vol(V{a}, Y{a});
                    
                    aas_log(aap,false,sprintf('Zero thresholded image %s sums up to %d vox', ...
                        V{a}.fname, sum(Y{a}(:))))
                end                
            case 'exclusive'
                %% C) Exclusive thresholding of masks:
                % Any particular voxel has greatest chance of being...
                maxeY = max(spm_read_vols(cell2mat(V)),[],4);
                
                for a = 1:size(segimg,1)
                    % We want to check where the segmentation has the
                    % greatest values
                    Y{a} = (Y{a} == maxeY) & (Y{a} > 0.01);
                    
                    V{a}.fname = spm_file(segimg(a,:),'prefix','E_r');
                    outstream = strvcat(outstream, V{a}.fname); % Save to stream...
                    spm_write_vol(V{a}, Y{a});
                    
                    aas_log(aap,false,sprintf('Exclusive thresholded image %s sums up to %d vox', ...
                        V{a}.fname, sum(Y{a}(:))))
                end
            otherwise
                %% B) Specific thresholding of each mask
                thr = aas_getsetting(aap,'threshold');
                if numel(thr) == 1, thr(1:3) = thr; end
                for a = 1:size(segimg,1)
                    maxY = max(Y{a}(:));
                    Y{a} = Y{a} > thr(a);
                    
                    if numel(Y{a} > 0) == 0 % check for bad thesholding
                        aas_log(aap, true, sprintf('ERROR: No voxels above the threshold mask (%f) [max: %f]', thr(a), maxY))
                    end
                    
                    V{a}.fname = spm_file(segimg(a,:),'prefix','S_r');
                    outstream = strvcat(outstream, V{a}.fname); % Save to stream...
                    spm_write_vol(V{a}, Y{a});
                    
                    aas_log(aap,false,sprintf('Strict thresholded image %s sums up to %d vox', ...
                        V{a}.fname, sum(Y{a}(:))))
                end                       
        end
        
        %% DIFFERENT STREAMS FOR DIFFERENT
        aap = aas_desc_outputs(aap,subj,'segmasks',outstream);
end
