% AA module - finetune coregistration of structural to mean EPI
% Should be done some time after normal (or extended) coregistration
% Coregistration of structural to mean EPI output by realignment in 3 steps
% 1) Copy the files to working copies...
% 2) Bias correct the mean EPI and structural
% 3) Rescale all the values remaining in the volumes from 0 to 255
% END) Do weighted coregistration, using the mean EPI mask
%


function [aap,resp]=aamod_coreg_finetune(aap,task,subj)

resp='';

switch task
    case 'doit'
        
        %% VARIOUS DEFAULTS
        
        % Coregistration
        global defaults
        flags = defaults.coreg;
        
        % Normalisation
        defs = aap.spm.defaults.normalise;
        defs.estimate.weight = '';
        
        % ...only write out attenuation corrected image
        estopts.regtype='';    % turn off affine:
        writeopts.biascor = 1;
        writeopts.GM  = [0 0 0];
        writeopts.WM  = [0 0 0];
        writeopts.CSF = [0 0 0];
        writeopts.cleanup = [0];
        
        % Realignment
        defs = aap.spm.defaults.realign;
        
        % ...flags to pass to routine to create resliced images
        % (spm_reslice)
        resFlags = struct(...
            'interp', defs.write.interp,...       % interpolation type
            'wrap', defs.write.wrap,...           % wrapping info (ignore...)
            'mask', defs.write.mask,...           % masking (see spm_reslice)
            'which', 1,...     % what images to reslice
            'mean', 0);           % write mean image
        
        %% 1) Get structural, mean EPI and masks...
        % Check local structural directory exists
        
        Simg = aas_getfiles_bystream(aap,subj,'structural');
        for b = 1:size(Simg, 1)
            % We want BETed image, ideally
            if strfind(Simg(b,:), 'bet_')
                Simg = deblank(Simg(b,:));
                break
            end
        end
        
        % Look for mean functional
        mEPIimg = aas_getfiles_bystream(aap,subj,1,'meanepi');
        if size(mEPIimg,1) > 1
            aas_log(aap, false, 'Found more than 1 mean functional images, using first.');
        end
        mEPIimg = deblank(mEPIimg(1,:));
        
        % Look for BET brain mask
        bSimg = aas_getfiles_bystream(aap,subj,'BETmask');
        for b = 1:size(bSimg, 1)
            if strfind(bSimg(b,:), 'brain_mask')
                bSimg = deblank(bSimg(b,:));
                break
            end
        end
        
        % Look for functional BET brain mask
        bmEPIimg = aas_getfiles_bystream(aap,subj,'epiBETmask');
        for b = 1:size(bmEPIimg, 1)
            if strfind(bmEPIimg(b,:), 'brain_mask')
                bmEPIimg = deblank(bmEPIimg(b,:));
                break
            end
        end
        
        if aap.tasklist.currenttask.settings.bias
            %% 2a) Bias correct images
            
            fprintf('Bias correct structural and mean EPI\n')
            
            [Spth, Sfn, Sext] = fileparts(Simg);
            [mEPIpth, mEPIfn, mEPIext] = fileparts(mEPIimg);
            
            out = spm_preproc(Simg, estopts);
            [sn,isn]   = spm_prep2sn(out);
            spm_preproc_write(sn, writeopts);
            
            out = spm_preproc(mEPIimg, estopts);
            [sn,isn]   = spm_prep2sn(out);
            spm_preproc_write(sn, writeopts);
            
            c_Simg = fullfile(Spth, ['m' Sfn Sext]);
            c_mEPIimg = fullfile(mEPIpth, ['m' mEPIfn mEPIext]);
        else
            %% 2b) Copy files to do preprocessing and coregistration on...
            
            fprintf('Copy images before processing them\n')
            
            [Spth, Sfn, Sext] = fileparts(Simg);
            [mEPIpth, mEPIfn, mEPIext] = fileparts(mEPIimg);
            
            c_Simg = fullfile(Spth, ['c_' Sfn Sext]);
            c_mEPIimg = fullfile(mEPIpth, ['c_' mEPIfn mEPIext]);
            
            unix(['cp ' Simg ' ' c_Simg])
            unix(['cp ' mEPIimg ' ' c_mEPIimg])
        end
        
        %% 3) RESCALE IMAGES...
        
        fprintf('Rescale the images to be coregistered\n')
        
        rescale4coreg(c_Simg)
        rescale4coreg(c_mEPIimg)
        
        %% 4)  Mean Functional to Structural (weighting with mEPI mask)
        
        % Add weighting...
        flags.estimate.wgt = spm_vol(bmEPIimg);
        
        % Estimate parameters
        x = spm_coreg_weighted(spm_vol(c_mEPIimg), ...
            spm_vol(c_Simg), ...
            flags.estimate);
        Mf = inv(spm_matrix(x));
        
        % Set the new space for the structural
        MM = spm_get_space(Simg);
        spm_get_space(Simg, Mf*MM);
        
        fprintf(['\tThe realignment parameters are the following\n' ...
            'x: %0.3f   y: %0.3f   z: %0.3f   p: %0.3f   r: %0.3f   j: %0.3f'], ...
            x(1), x(2), x(3), x(4), x(5), x(6))
        
        %% Some diagnostic images
        spm_check_registration(strvcat( ...
            Simg, ... % Get structural
            mEPIimg)); % Get mean EPI across sessions
        
        % Save graphical output to common diagnostics directory
        if ~exist(fullfile(aap.acq_details.root, 'diagnostics'), 'dir')
            mkdir(fullfile(aap.acq_details.root, 'diagnostics'))
        end
        mriname = strtok(aap.acq_details.subjects(subj).mriname, '/');
        figure(spm_figure('FindWin'));
        set(gcf,'PaperPositionMode','auto')
        print('-djpeg','-r75',fullfile(aap.acq_details.root, 'diagnostics', ...
            [mfilename '__' mriname '.jpeg']));
        
        %% Diagnostic VIDEO of coregistration
        
        if aap.tasklist.currenttask.settings.diagnostic
            
            Ydims = {'X', 'Y', 'Z'};
            
            % Get mean EPI
            Y = spm_read_vols(spm_vol(mEPIimg));
            
            % Get resliced structural
            warning off
            spm_reslice(strvcat(mEPIimg, Simg), resFlags);
            sY = spm_read_vols(spm_vol(fullfile(Spth, ['r' Sfn Sext])));
            warning on
            
            EPIlims = [min(Y(:)) max(Y(:))];
            
            for d = 1:length(Ydims)
                movieFilename = fullfile(aap.acq_details.root, 'diagnostics', ...
                    [mfilename '__' mriname '_' Ydims{d} '.avi']);
                % Create movie file by defining aviObject
                try delete(movieFilename); catch; end
                aviObject = avifile(movieFilename,'compression','none');
                
                try close(2); catch; end
                figure(2)
                set(2, 'Position', [0 0 1000 800])
                windowSize = get(2,'Position');
                
                for n = 1:size(sY,d)
                    % Get outline of structural slice
                    h = subplot(1,1,1);
                    if d == 1
                        sOutline = edge(rot90(squeeze(sY(n,:,:))),'canny');
                    elseif d == 2
                        sOutline = edge(rot90(squeeze(sY(:,n,:))),'canny');
                    elseif d == 3
                        sOutline = edge(rot90(squeeze(sY(:,:,n))),'canny');
                    end
                    
                    % Get image of EPI slice
                    if d == 1
                        sImage = rot90(squeeze(Y(n,:,:)));
                    elseif d == 2
                        sImage = rot90(squeeze(Y(:,n,:)));
                    elseif d == 3
                        sImage = rot90(squeeze(Y(:,:,n)));
                    end
                    
                    % Overlay structural outline on EPI image
                    sImage(logical(sOutline)) = EPIlims(2) * 2;
                    imagesc(sImage)
                    
                    caxis(EPIlims)
                    axis equal off
                    zoomSubplot(h, 1.2)
                    
                    % Capture frame and store in aviObject
                    pause(0.01)
                    aviObject = addframe(aviObject,getframe(2,windowSize));
                end
                
                aviObject = close(aviObject);
            end
            try close(2); catch; end
        end
        
        %% Describe the outputs
        
        aap = aas_desc_outputs(aap,subj,'structural',Simg);
        
    case 'checkrequirements'
        aas_log(aap,0,'No need to trim or skull strip structural\n' );
end
