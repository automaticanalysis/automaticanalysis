% AA module
% High pass filter a timeseries of data

function [aap,resp]=aamod_highpass(aap,task,subj,sess)

resp='';

switch task
    case 'domain'
        resp='session';  % this module needs to be run once per subject
        
    case 'description'
        resp='SPM highpass filter';
        
    case 'summary'
        subjpath=aas_getsubjpath(subj);
        resp=sprintf('\n\tHighpass %s\n',subjpath);
        
    case 'report'
        
    case 'doit'
        
        EPIimg = aas_getfiles_bystream(aap,subj,sess,'epi');
        
        %% retrieve TR from DICOM header & set up the HiPass filter
        % if TR is manually specified (not recommended as source of error)
        if (isfield(aap.tasklist.currenttask.settings,'TR'))
            K.RT =aap.tasklist.currenttask.settings.TR;
        else
            % Get TR from DICOM header
            DICOMHEADERS=load(aas_getfiles_bystream(aap,subj,sess,'epi_dicom_header'));
            K.RT = DICOMHEADERS.DICOMHEADERS{1}.volumeTR;
        end
        
        % High pass filter or detrend data
        % Let's first set up the parameters...
        K.row = 1:size(EPIimg, 1);
        K.HParam = aap.tasklist.currenttask.settings.HParam; % cut-off period in seconds
        
        if K.RT * length(K.row) > K.HParam
            fprintf('\nWill do high pass filtering of time series with a %d second cut-off', K.HParam)
        else
            fprintf('\nWill do linear detrending across time series')
        end
        
        V = spm_vol(deblank(EPIimg(1,:))); % A typical volume...
        
        %% If the dataset is too large, we process it by chunks...
        fprintf('\n\tProcessing data (%d scans)', size(EPIimg,1))
        
        taskComplete = 0;
        chunkDim = 1;
        
        while taskComplete == 0
            fprintf('\nTrying with %d chunks', chunkDim)
            
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
                            fprintf('\n\t...chunk %d %d %d', x, y, z)
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
                            end
                            
                            if K.RT * length(K.row) > K.HParam
                                % Create the frequencies to be removed and apply them...
                                % Important: first dimension must be time dimension!
                                
                                EPIdata = spm_filter(K, EPIdata);
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
                            end
                            
                            % Now save the data back...
                            for e = 1:size(EPIimg,1)
                                V = spm_vol(deblank(EPIimg(e,:)));
                                Y = spm_read_vols(V);
                                Y(Xind,Yind,Zind) = EPIdata(e,:,:,:);
                                spm_write_vol(V,Y);
                            end
                        end
                    end
                end
                % If we get here, then we completed the task...
                taskComplete = 1;
            catch tSNR_error
                %disp(tSNR_error)
                
                if x > 1 || y > 1 || z > 1
                    aas_log(aap, true, 'The script broke between chunks, you should probably delete the subject folder for this module and try again...')
                end
                
                if chunkDim > 3
                    aas_log(aap, true, 'Error is probably not due to MEMORY')
                end
                
                chunkDim = chunkDim + 1;
            end
        end
        
        %% DESCRIBE THE OUTPUTS
        
        aap=aas_desc_outputs(aap,subj, sess, 'epi', EPIimg);
        
end
