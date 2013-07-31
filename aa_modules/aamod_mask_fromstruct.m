% AA module - Mask images from a segmented structural
% 1) Take segmentation images origination from unified segmentation (and
% normalisation)
% 2) Threshold them at 3 different levels (zero level, strict [e.g. 99%],
% or an exclusive [highest tissue probability wins])

function [aap,resp] = aamod_mask_fromstruct(aap,task,subj)

resp = '';

switch task
    case 'doit'
                
        Simg = aas_getfiles_bystream(aap,subj,'structural');
        SEGimg = aas_getfiles_bystream(aap,subj,'segmentation');
        mEPIimg = aas_getfiles_bystream(aap,subj,1,'meanepi');
        
        % Whether the images in SEGimg are native (1) or warped (0)
        NWlogical = zeros(size(SEGimg,1), 1);
        
        % Cheap and cheerful way of ensuring only one file is considered!
        if size(Simg,1) > 1
            for a = 1:size(Simg,1)
                % Not warped or betted!
                if ~strcmp(Simg(a,1), 'w') && ~strcmp(Simg(a,1), 'b')
                    Simg = deblank(Simg(a,:));
                    break
                end
            end
            fprintf('\tSeveral structurals found, considering: %s\n', Simg)
        end
        if size(mEPIimg,1) > 1
            % Not warped!
            for a = 1:size(mEPIimg,1)
                if ~strcmp(mEPIimg(a,1), 'w')
                    mEPIimg = deblank(mEPIimg(a,:));
                    break
                end
            end
            mEPIimg = mEPIimg(1,:);
            fprintf('\tSeveral mean EPIs found, considering: %s\n', mEPIimg)
        end
        
        %% 1) RESLICE
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
        
        for a = 1:size(SEGimg,1)
            [junk, fn] = fileparts(SEGimg(a,:));
            if strcmp(fn(1), 'w')
                % Reslice warped segmentations to  warped structural image
                spm_reslice(strvcat(Simg,SEGimg(a,:)))
            else
                % Reslice native segmentations to the mean EPI
                spm_reslice(strvcat(mEPIimg,SEGimg(a,:)))
                NWlogical(a) = 1; % Set this as a native image...
            end
        end
        
        %% 2) THRESHOLD
        
        %% IN PROGRESS
        Zoutstream = '';
        Soutstream = '';
        Eoutstream = '';
        
        NeY = 0;
        WeY = 0;
        maxWeY = 0;
        maxNeY = 0;
        
        V = cell(1, size(SEGimg,1));
        Y = cell(1, size(SEGimg,1));
        
        %% A) Zero thresholding is easy...
        for a = 1:size(SEGimg,1)
            [pth, fn, ext] = fileparts(SEGimg(a,:));
            % Load the correct resliced file!
            V{a} = spm_vol(fullfile(pth, ['r' fn ext]));
            Y{a} = spm_read_vols(V{a});
            
            zY = Y{a} > 0;
            V{a}.fname = fullfile(pth, ['Z_r' fn ext]);
            Zoutstream = strvcat(Zoutstream, V{a}.fname); % Save to stream...
            spm_write_vol(V{a}, zY);
            
            fprintf('Zero thresholded image %s sums up to %d vox\n', ...
                fn, sum(zY(:)))
        end
        fprintf('\n')
        
        %% B) Specific thresholding of each mask
        for a = 1:size(SEGimg,1)
            [pth, fn, ext] = fileparts(SEGimg(a,:));
            
            if strcmp(fn(1), 'w')
                % Warped threshold
                sY = Y{a} > aap.tasklist.currenttask.settings.Wthreshold;
            else
                % Native threshold
                sY = Y{a} > aap.tasklist.currenttask.settings.Nthreshold;
            end
            V{a}.fname = fullfile(pth, ['S_r' fn ext]);
            Soutstream = strvcat(Soutstream, V{a}.fname); % Save to stream...
            spm_write_vol(V{a}, sY);
            
            fprintf('Strict thresholded image %s sums up to %d vox\n', ...
                fn, sum(sY(:)))
            
            %% C) Exclusive thresholding of masks:
            % Any particular voxel has greatest chance of being...
            
            % We want to check where the new segmentation has greater
            % values than the last time
            if ~isempty(strfind(fn(1), 'w'))
                if size(WeY,2) == 1;
                    WeY = zeros(size(Y{a}));
                end
                % Warped threshold
                WeY(Y{a} > maxWeY) = a;
                maxWeY = max(Y{a}, maxWeY);
            else
                if size(NeY,2) == 1;
                    NeY = zeros(size(Y{a}));
                end
                % Native threshold
                NeY(Y{a} > maxNeY) = a;
                maxNeY = max(Y{a}, maxNeY);
            end
        end
        fprintf('\n')
        
        % Last stage of exclusive thresholding
        for a = 1:size(SEGimg,1)
            [pth, fn, ext] = fileparts(SEGimg(a,:));
            
            if ~isempty(strfind(fn(1), 'w'))
                eY = WeY == a;
            else
                eY = NeY == a;
            end
            eY = double(eY);
            
            V{a}.fname = fullfile(pth, ['E_r' fn ext]);
            Eoutstream = strvcat(Eoutstream, V{a}.fname); % Save to stream...
            spm_write_vol(V{a}, eY);
            
            fprintf('Exclusive thresholded image %s sums up to %d vox\n', ...
                fn, sum(eY(:)))
        end
        
        % DIFFERENT STREAMS FOR DIFFERENT
        aap = aas_desc_outputs(aap,subj,'segmasksZero',Zoutstream);
        aap = aas_desc_outputs(aap,subj,'segmasksStrict',Soutstream);
        aap = aas_desc_outputs(aap,subj,'segmasksExclusive',Eoutstream);
        
end