% AA module
% [aap,resp]=aamod_tSNR_EPI(aap,task,subj,sess)
% Calculate the temporal SNR of a series of images within a particular ROI
% Also plots the timecourses...
% This is more accurate (but slower!) than the SNR calculation that is done
% during epi conversion, since it removes signal drift!

function [aap,resp]=aamod_tSNR_EPI(aap,task,subj,sess)

resp='';

switch task
    case 'summary'
        sesspath=aas_getsesspath(aap,subj,sess);
        resp=sprintf('Estimate tSNR of EPI series of %s\n',sesspath);
        
    case 'report'
        
    case 'doit'
        
        EPIimg = aas_getfiles_bystream(aap,subj,sess,'epi');
        if aap.options.NIFTI4D
            EPIimg = spm_select('ExtFPList',fileparts(EPIimg),basename(EPIimg));
        end
        if aas_stream_has_contents(aap,'rois')
            ROIimg = aas_getfiles_bystream(aap,subj,'rois');
        else
            ROIimg = [];
        end
        
        %% retrieve TR from DICOM header & set up the HiPass filter
        K.RT = [];
        % if TR is manually specified (not recommended as source of error)
        if isfield(aap.tasklist.currenttask.settings,'TR') && ~isempty(aap.tasklist.currenttask.settings.TR)
            K.RT = aap.tasklist.currenttask.settings.TR;
        else
            % Get TR from DICOM header
            if aas_stream_has_contents(aap,subj,sess,'epi_dicom_header');
                DICOMHEADERS = load(aas_getfiles_bystream(aap,subj,sess,'epi_dicom_header'));
                K.RT = DICOMHEADERS.DICOMHEADERS{1}.volumeTR;
            end
        end
        % High pass filter or detrend data
        % Let's first set up the parameters...
        K.row = 1:size(EPIimg,1);
        K.HParam = aap.tasklist.currenttask.settings.HParam; % cut-off period in seconds
        % If total length > filter cut-off
        
        doFilt = ~isempty(K.RT) && ~isempty(K.HParam);
        
        if doFilt
            if K.RT * length(K.row) > K.HParam
                aas_log(aap,false,sprintf('Will do high pass filtering of time series with a %d second cut-off', K.HParam));
            else
                aas_log(aap,false,'Will do linear detrending across time series')
            end
        else
            aas_log(aap,false,'No high pass filtering will be applied on the time series')
        end
        
        %% Get started with the processing
        % Voxel based measures
        V = spm_vol(deblank(EPIimg(1,:))); % A typical volume...
        EPIsignal = zeros(V.dim(1), V.dim(2), V.dim(3));
        EPInoise = zeros(V.dim(1), V.dim(2), V.dim(3));
        
        if ~isempty(ROIimg)
            % ROIvol{r}based measures
            EPIsnROI= cell(size(ROIimg,1),1);
            EPIsnHist = cell(size(ROIimg,1),1);
            mROI = cell(size(ROIimg,1),1);
            sROI = cell(size(ROIimg,1),1);
            SNmROI= cell(size(ROIimg,1),1);
            ROIname = cell(size(ROIimg,1),1);
            ROIvol = cell(size(ROIimg,1),1);
        end
        
        aas_log(aap,false,sprintf('Working on session %s', aap.acq_details.sessions(sess).name))
        
        aas_log(aap,false,'\tLoading ROIs')
        for r = 1:size(ROIimg,1)
            % Now load each of the ROIs we wish to examine (usually using the grey matter)
            rV = spm_vol(ROIimg(r,:));
            ROIvol{r} = spm_read_vols(rV);
            ROIvol{r} = round(ROIvol{r});
            % Mean ROI value is a vector depending on number of scans
            mROI{r} = zeros(size(EPIimg,1),1);
            
            if any(size(ROIvol{r})~=size(EPIsignal))
                aas_log(aap, true, ['The dimensions of the EPI data and the ROI do not match\n' ...
                    'This is likely because you are using normalised EPIs, whereas you should use native ones'])
            end
        end
        
        %% If the dataset is too large, we process it by chunks...
        aas_log(aap,false,sprintf('\tProcessing data (%d scans)', size(EPIimg,1)))
        
        taskComplete = 0;
        chunkDim = 1;
        
        while taskComplete == 0
            aas_log(aap,false,sprintf('Trying with %d chunks', chunkDim))
            
            try
                chunkX = 0;
                chunkY = 0;
                chunkZ = 0;
                for c = 1:chunkDim
                    chunkX = [chunkX floor(V.dim(1) * c / chunkDim)];
                    chunkY = [chunkY floor(V.dim(2) * c / chunkDim)];
                    chunkZ = [chunkZ floor(V.dim(3) * c / chunkDim)];
                end
                
                % Chunking...
                for x = 1:length(chunkX) - 1
                    for y = 1:length(chunkY) - 1
                        for z = 1:length(chunkZ) - 1
                            aas_log(aap,false,sprintf('\t...chunk %d %d %d', x, y, z))
                            Xind = chunkX(x) + 1 : chunkX(x+1);
                            Yind = chunkY(y) + 1 : chunkY(y+1);
                            Zind = chunkZ(z) + 1 : chunkZ(z+1);
                            
                            EPIdata = zeros(size(EPIimg,1), length(Xind), ...
                                length(Yind), ...
                                length(Zind));
                            
                            % Load each image into 4-D matrix
                            for e = 1:size(EPIimg,1)
                                V = spm_vol(deblank(EPIimg(e,:)));
                                Y = spm_read_vols(V);
                                EPIdata(e,:,:,:) = Y(Xind,Yind,Zind);
                                
                                %% We can do ROIvol{r} processing here...
                                if x == 1 && y == 1 && z == 1
                                    for r = 1:size(ROIimg,1)
                                        tmp = Y(ROIvol{r}>0);
                                        tmp = tmp(tmp>0); % We don't want to include zero values...
                                        mROI{r}(e,1) = mean(tmp(:)); % Mean per time point
                                        sROI{r}(e,1) = std(tmp(:))./sqrt(length(tmp(:))); % Standard error per time point
                                    end
                                end
                            end
                            
                            if doFilt
                                if K.RT * length(K.row) > K.HParam
                                    % Create the frequencies to be removed and apply them...
                                    % Important: first dimension must be time dimension!
                                    EPIdata = reshape(spm_filter(K, reshape(EPIdata,size(EPIdata,1),[])),size(EPIdata));
                                    % And the ROI data...
                                    if x == 1 && y == 1 && z == 1
                                        for r = 1:size(ROIimg,1)
                                            mROI{r}(:,1) = spm_filter(K, mROI{r});
                                        end
                                    end
                                else
                                    % Use linear detrending instead (might be slower due to loops)
                                    for a = 1:length(Xind)
                                        for b = 1:length(Yind)
                                            vRow = squeeze(EPIdata(:,a,b,:));
                                            mRow = repmat(mean(vRow,1), [size(EPIimg,1) 1]);
                                            vRow = detrend(vRow);
                                            % Add mean back after detrending!
                                            EPIdata(:,a,b,:) = vRow + mRow;
                                        end
                                    end
                                    % And the ROI data...
                                    if x == 1 && y == 1 && z == 1
                                        for r = 1:size(ROIimg,1)
                                            mROI{r}(:,1) = detrend(mROI{r}) + mean(mROI{r});
                                        end
                                    end
                                end
                            end
                            
                            % Calcultate signal as mean of the data across volumes
                            EPIsignal(Xind,Yind,Zind) = squeeze(mean(EPIdata, 1));
                            % Calculate noise as standard deviation across volumes
                            EPInoise(Xind,Yind,Zind) = squeeze(std(EPIdata, [], 1));
                        end
                    end
                end
                % If we get here, then we completed the task...
                taskComplete = 1;
            catch tSNR_error
                msg = sprintf('Error: %s\n',tSNR_error.message);
                for e = 1:numel(tSNR_error.stack)
                    msg = [msg sprintf('in %s (line %d)\n', ...
                        tSNR_error.stack(e).file, tSNR_error.stack(e).line)];
                end
                
                if chunkDim > 3
                    msg = [msg 'Error is probably not due to MEMORY\n'];
                end
                aas_log(aap,false,msg,aap.gui_controls.colours.error)
                chunkDim = chunkDim + 1;
            end
        end
        
        aas_log(aap,false,'\tCalculating & saving the tSNR image')
        % Calculate SNR as ratio of the two...
        EPIsnr = EPIsignal ./ EPInoise;
        EPIsnr(isnan(EPIsnr)|isinf(EPIsnr)) = 0;
        
        % Avoids extreeme outliers!
        EPImaxSNR = squeeze(max(EPIsnr(EPIsnr<1000)));
        
        % Save the SNR image!
        sV = V;
        sV.n = [1 1]; % make sure it is 3D
        sV.fname = fullfile(aas_getsesspath(aap,subj,sess), ...
            ['tSNR_' aap.acq_details.sessions(sess).name '.nii']);
        spm_write_vol(sV, EPIsnr);
        
        aas_log(aap,false,'\tFinalising ROI data')
        
        for r = 1:size(ROIimg,1)
            % Now get the voxels specific to each ROI
            EPIsnROI{r} = EPIsnr(ROIvol{r}>0);
            EPIsnROI{r} = EPIsnROI{r}(EPIsnROI{r}>0); % We don't want to include zero values...
            % Also get a whole ROI signal and noise estimate
            SNmROI{r} = mean(mROI{r}) ./ std(mROI{r});
        end
        
        %% DIAGNOSTIC IMAGE
        % Save graphical output to common diagnostics directory
        if ~exist(fullfile(aap.acq_details.root, 'diagnostics'), 'dir')
            mkdir(fullfile(aap.acq_details.root, 'diagnostics'))
        end
        subjname = aap.acq_details.subjects(subj).subjname;
        for r = 1:size(ROIimg,1)
            [junk, ROIname{r}] = fileparts(ROIimg(r,:));
            
            %% Diagnostic VIDEO of masks
            if aap.tasklist.currenttask.settings.diagnostic
                
                movieFilename = fullfile(aap.acq_details.root, 'diagnostics', ...
                    [mfilename '__' subjname '_' ROIname{r} '.avi']);
                % Create movie file by defining aviObject
                try delete(movieFilename); catch; end
                aviObject = avifile(movieFilename,'compression','none');
                windowSize = [0 0 800 800];
                
                % Plot the ROIvol{r}overlaid onto mean data
                try close(2); catch; end
                figure(2)
                set(2, 'Position', windowSize)
                
                for x = 1:size(ROIvol{r},1)
                    h = subplot(1,1,1);
                    % Edge of the ROI
                    sOutline = rot90(edge(squeeze(ROIvol{r}(x,:,:)), 'canny'));
                    
                    sImage = rot90(squeeze(EPIsnr(x,:,:)));
                    sImage(sOutline) = EPImaxSNR * 2;
                    
                    imagesc(sImage);
                    caxis([0 EPImaxSNR])
                    colorbar
                    axis equal off
                    title(sprintf('SNR of our scan, overlaid with ROI: %s',ROIname{r}))
                    zoomSubplot(h, 1.2)
                    
                    pause(0.01)
                    % Capture frame and store in aviObject
                    aviObject = addframe(aviObject,getframe(2,windowSize));
                end
                
                junk = close(aviObject);
            end
        end
        
        %% tSNR results figure!
        if ~isempty(ROIimg)
            aas_log(aap,false,'Displaying the results of the tSNR analysis')
            colorsB = {'r' 'g' 'b' 'c' 'm' 'y' 'w'};
            
            % We need to make a string for legend...
            legStr = cell(1,size(ROIimg,1));
            for r = 1:size(ROIimg,1)
                legStr{r} = sprintf('%s; mn:%.2f; SD:%.2f; med:%.0f; ROI:%.2f (%.0fv)', ...
                    ROIname{num2str(r)}, ...
                    mean(EPIsnROI{num2str(r)}), ...
                    std(EPIsnROI{num2str(r)}), ...
                    median(EPIsnROI{num2str(r)}), ...
                    SNmROI{num2str(r)}, ...
                    sum(EPIsnROI{num2str(r)}(:)));
            end
            
            try close(2); catch; end
            
            figure(2)
            set(2, 'Position', [0 0 1200 700])
            maxI = 0;
            windI = 0;
            maxV = 0;
            hold on
            
            for r = 1:size(ROIimg,1)
                % What range do the SNR values take?
                maxI = max(max(EPIsnROI{r}), maxI);
                % What window do we wish to present?
                windI = max(median(EPIsnROI{r}) + std(EPIsnROI{r}) * 3, windI);
            end
            vals = 0:windI/100:ceil(maxI);
            for r = 1:size(ROIimg,1)
                % Now make a histogram and "normalise" it
                EPIsnHist{r} = hist(EPIsnROI{r}, vals);
                EPIsnHist{r} = EPIsnHist{r}./sum(EPIsnHist{r});
                % And decide what is the greatest prop value
                maxV = max(max(EPIsnHist{r}), maxV);
                
                % Make bars semi-transparent for overlaying
                B = bar(vals, EPIsnHist{r}, 1, colorsB{r});
                ch = get(B,'child');
                set(ch, 'faceA', 0.3, 'edgeA', 0.2);
            end
            vals = 0:windI/100:ceil(windI);
            % Set the axis to a good value!
            axis([vals(1), vals(end), 0, maxV*1.05])
            xlabel('SNR')
            ylabel('Proportion of voxels')
            set(gca,'XTick', 0:ceil(maxI./50):maxI)
            title(sprintf('\nSNR for session %s, using %.0f scans and TR %.3f', ...
                regexprep(aap.acq_details.sessions(sess).name, '[^a-zA-Z0-9]', ''), ...
                size(EPIimg,1), ...
                K.RT))
            legend(legStr);
            
            set(gcf,'PaperPositionMode','auto','Renderer','zbuffer');
            print('-djpeg','-r75',fullfile(aap.acq_details.root, 'diagnostics', ...
                [mfilename '__' subjname '_tSNR.jpeg']));
            
            %% Time-course results figure!
            
            aas_log(aap,false,'Displaying the results of the timecourse analysis')
            
            try close(2); catch; end
            
            figure(2)
            set(2, 'Position', [0 0 1200 700])
            
            for r = 1:size(ROIimg,1)
                subplot(size(ROIimg,1),1,r)
                hold on
                
                % We need to make a string for legend...
                legStr = sprintf('%s (%.0fv)', ...
                    ROIname{num2str(r)}, ...
                    sum(EPIsnROI{num2str(r)}(:)));
                
                % Plot main results (errorbars displayed differently now...)
                plot(mROI{r}, ['.' colorsB{r}])
                plot(mROI{r} - sROI{r}, '--k')
                plot(mROI{r} + sROI{r}, '--k')
                
                xlabel('Scan')
                ylabel('Mean signal')
                set(gca,'XTick', 0:ceil(size(EPIimg,1)./25):size(EPIimg,1))
                xlim([0 size(EPIimg,1)])
                ylim([mean(mROI{r} - 2*mean(sROI{r})) mean(mROI{r} + 2*mean(sROI{r}))])
                title(sprintf('\nTimecourse for session %s, using %.0f scans and TR %.3f', ...
                    regexprep(aap.acq_details.sessions(sess).name, '[^a-zA-Z0-9]', ''), ...
                    size(EPIimg,1), ...
                    K.RT))
                legend(legStr);
            end
            
            set(gcf,'PaperPositionMode','auto','Renderer','zbuffer');
            print('-djpeg','-r75',fullfile(aap.acq_details.root, 'diagnostics', ...
                [mfilename '__' subjname '_timecourse.jpeg']));
        end
        %% DESCRIBE OUTPUTS
        
        aap=aas_desc_outputs(aap,'session',[subj sess],'tSNR',sV.fname);
end