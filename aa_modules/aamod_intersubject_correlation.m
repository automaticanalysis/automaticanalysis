function [aap,resp] = aamod_intersubject_correlation(aap, task)
%
% generate intersubject correlation maps from epi
%
% Change History
%
% spring 2023 [MSJ] -- new
%

resp='';

switch task
	
    case 'report'
        
    case 'doit'		
        
        % put some xml settings into local variables for better code readability

        subtract_global_signal = aap.tasklist.currenttask.settings.subtract_global_signal;
        convert_r_to_z = aap.tasklist.currenttask.settings.convert_r_to_z;       
        verbose = aap.tasklist.currenttask.settings.verbose;
        save_individual_maps = aap.tasklist.currenttask.settings.save_individual_maps;
        
        rmap_threshold = aap.tasklist.currenttask.settings.rmap_threshold;
        show_only_positive_values = aap.tasklist.currenttask.settings.show_only_positive_values; 
        show_only_negative_values = aap.tasklist.currenttask.settings.show_only_negative_values; 
        
        maximum_number_of_sessions = aap.tasklist.currenttask.settings.maximum_number_of_sessions;

        if (save_individual_maps == true)
            mkdir(aas_getstudypath(aap),fullfile('pairwise_maps'));
            individual_map_fnames = {}; % cell array of fnames for desc
        end
        
        generate_summary_correlation_matrix = aap.tasklist.currenttask.settings.generate_summary_correlation_matrix;
        
        % ****************     
        % DATA
        % ****************
        
        % assemble a big list of nii files (all subjects, all sessions)
        
        fnames = {};
        FID = {};

        [ ~,subindex_list ] = aas_getN_bydomain(aap,'subject');

        for subject_index = subindex_list

            [ ~,all_sess ] = aas_getN_bydomain(aap, 'session', subject_index);
            session_list = intersect(all_sess, aap.acq_details.selected_sessions);           
            
            % clip to number of sessions requested if necessary
            
           if (length(session_list) > maximum_number_of_sessions) 
                session_list = session_list(1:maximum_number_of_sessions);
            end

            for session_index = session_list

                temp = aas_getfiles_bystream(aap, subject_index, session_index, 'epi');              
                fnames{end+1} = temp;
                FID{end+1} = sprintf('SUB%02dSESS%02d',subject_index, session_index);

            end

        end 
        
        % ****************     
        % MASKING
        % ****************
                       
        % three options are provided for specifying a mask:
        %
        %   1) a full path to a file in <explicit_brain_mask_fname>
        %   2) specified using the optional renamable input stream "brainmask"
        %   3) on-the-fly implicit mask using threshold setting
        %
        % if more than one mask is defined, the intersection is applied
                 
        % the bare minimum mask simply excludes zero and NaN voxels 
        % we're going to assume we can determine this from the first epi
                     
        header = spm_vol(fnames{1});
        Y = spm_read_vols(header);             % nrow x ncol x nslice x nframes       
        mask = all( (Y~=0 & ~isnan(Y)), 4);    % ncow x ncol x nslice (logical)
                       
        % we intersect minimal mask with implicit, instream and/or explicit masks 
        % masks (whichever are defined) to generate final mask
                
        % 1) optional implicit mask
        
        % assume we can build this using first epi (currently in "Y")
        
        if (~isempty(aap.tasklist.currenttask.settings.implicit_mask_threshold))
            implicit_mask_threshold = aap.tasklist.currenttask.settings.implicit_mask_threshold;          
            mean_over_t = mean(abs(Y),4,'omitnan'); % nrow x ncol x nslice
            global_mean = mean(mean_over_t(:)); % scalar
            implicit_mask = mean_over_t > implicit_mask_threshold*global_mean;
            mask = mask & implicit_mask;
        end
 
        % we're done with (big variable) Y
        
        clear Y;
        
        % save mask now in case reslice is needed for other options

        mask_header = header(1);
        mask_fname = fullfile(aas_getstudypath(aap), 'mask.nii');
        mask_header.fname = mask_fname;
        mask_header.descrip = 'ISC mask';
        spm_write_vol(mask_header, mask);
        pause(1); % let filesystem catch up
        
        % but may need to save again if there's an intersection
        
        mask_has_changed = false; 
        
        % 2) optional instream mask

        mask_inputstream_struct = aap.tasklist.currenttask.inputstreams(1).stream{1};

        if (aas_stream_has_contents(aap, mask_inputstream_struct.CONTENT))

            instream_mask_fname = aas_getfiles_bystream(aap, mask_inputstream_struct.CONTENT);
            instream_mask_header = spm_vol(instream_mask_fname);

            if (~isequal(instream_mask_header.dim,mask_header.dim)|| ~isequal(instream_mask_header.mat,mask_header.mat))

               temp_fname = reslice_file_to_file(mask_fname, instream_mask_fname);
               instream_mask_header = spm_vol(temp_fname);        

            end

            instream_mask = spm_read_vols(instream_mask_header);
            mask = mask & instream_mask;
            mask_has_changed = true;
            
        end

        % 3) optional explicit mask
        
        if (~isempty(aap.tasklist.currenttask.settings.explicit_mask_fname))

            explicit_mask_fname = aap.tasklist.currenttask.settings.explicit_mask_fname;
            explicit_mask_header = spm_vol(explicit_mask_fname);

            if (~isequal(explicit_mask_header.dim, mask_header.dim)|| ~isequal(explicit_mask_header.mat,mask_header.mat))

                temp_fname = reslice_file_to_file(mask_fname, explicit_mask_fname);
                explicit_mask_header = spm_vol(temp_fname);        

            end

            explicit_mask = spm_read_vols(explicit_mask_header);
            mask = mask & explicit_mask;
            mask_has_changed = true; 

        end

        % if we updated the mask, need to save it again...
        
        if (mask_has_changed) 
            mask_fname = fullfile(aas_getstudypath(aap), 'mask.nii');
            mask_header.fname = mask_fname;
            spm_write_vol(mask_header, mask);      
        end
        
        % save jpgs of mask for QA
        
        save_three_ortho_jpgs(mask,fullfile(aas_getstudypath(aap),'BRAINMASK'));       
        
        % we use the mask header as a template for all subsequent file saves
        
        template_header = mask_header;
        
        % several of the options need an "unrolled" version of the mask
        % create it here, once
        
        [nrow,ncol,nslice] = size(mask);
        unrolled_mask = reshape(mask,nrow*ncol*nslice,1);

        % -----------------------------------------------------------------
        % loop over all files in fnames, generates pairwise corr maps 
        % -----------------------------------------------------------------
                                      
        if (generate_summary_correlation_matrix == true)
            summary_correlation_matrix_r = eye(numel(fnames));
            summary_correlation_matrix_p = zeros(numel(fnames));
            summary_descriptors = {};
        end

        file_count = 0;  % tally for averaging

        for findex1 = 1:numel(fnames)-1

            fname1 = fnames{findex1};  
            header = spm_vol(fname1);
            data1 = spm_read_vols(header);

            [nrow,ncol,nslice,nframe] = size(data1);
            
            if (subtract_global_signal == true)
                tempdata = reshape(data1,nrow*ncol*nslice,nframe);
                tempdata = tempdata - mean(tempdata(unrolled_mask==1,:),1);
                % subtraction may have created nonzero values
                % outside of mask; ergo, reapply mask now:
                tempdata(unrolled_mask~=1,:) = 0;
                data1 = reshape(tempdata,nrow,ncol,nslice,nframe);
            end

            if (generate_summary_correlation_matrix == true)
                % extract summary data *before* normalization
                tempdata1 = reshape(data1,nrow*ncol*nslice,nframe);
                data1_summary = mean(tempdata1(unrolled_mask==1,:),1);
            end
            
            data1 = normalize(data1,4);
    
            for findex2 = findex1+1:numel(fnames)

                if (verbose)
                    aas_log(aap,false,sprintf('INFO: processing %s v %s', FID{findex1}, FID{findex2}));
                end

                fname2 = fnames{findex2};        
                header = spm_vol(fname2);
                data2 = spm_read_vols(header); 
          
                if (subtract_global_signal == true)
                    tempdata = reshape(data2,nrow*ncol*nslice,nframe);
                    tempdata = tempdata - mean(tempdata(unrolled_mask==1,:),1);
                    tempdata(unrolled_mask~=1,:) = 0;
                    data2 = reshape(tempdata,nrow,ncol,nslice,nframe);
                end 

               if (generate_summary_correlation_matrix == true)
                    tempdata2 = reshape(data2,nrow*ncol*nslice,nframe);
                    data2_summary = mean(tempdata2(unrolled_mask==1,:),1); 
                    [ summary_r,summary_p ] = corr(data1_summary(:), data2_summary(:));
                    summary_correlation_matrix_r(findex1,findex2) = summary_r;
                    summary_correlation_matrix_p(findex1,findex2) = summary_p;
                    summary_descriptors{end+1} =  sprintf('%s_v_%s', FID{findex1}, FID{findex2});                 
                end
                
                data2 = normalize(data2,4);

                % for normalized data, r = dot product divided by vector length - 1

                RMAP = dot(data1,data2,4) / (size(data1,4)-1);

                % normalization will create NaN at voxels outside of brain
                % (because these voxels have constant (i.e. 0) timeseries)
                % -- reapply mask to zero these

                RMAP(mask==0) = 0;
                
                nancheck = sum(isnan(RMAP(:)));
                if (verbose && nancheck > 0)
                    aas_log(aap,false,sprintf('INFO: %d NaN survive after remasking.', nancheck));
                end                

                if (verbose)
                    aas_log(aap,false,sprintf('INFO: rmin: %6.3f rmax: %6.3f', min(RMAP(:)), max(RMAP(:))));
                end
 
                % threshold and clipping
                
                if (rmap_threshold > -inf)
                    RMAP(RMAP<rmap_threshold) = 0;
                end
                
                if (show_only_positive_values == true)
                    RMAP(RMAP<0) = 0;
                end
                
               if (show_only_negative_values == true)
                    RMAP(RMAP>0) = 0;
                end

                % add this map for later averaging; note we
                % sum Z transformed maps (then untransform at end)
 
                if (file_count == 0)
                    SUMMAP = atanh(RMAP);
                    file_count = 1;
                else
                    SUMMAP = SUMMAP + atanh(RMAP);
                    file_count = file_count + 1; 
                end

                % save the individual map?

                if (save_individual_maps == true)
                    
                    if (convert_r_to_z)
                        RMAP = atanh(RMAP);
                        template_header.descrip = 'ISC (z values)';
                    else
                        template_header.descrip = 'ISC (r values)';
                    end

                    temp = sprintf('%s_%s.nii', FID{findex1}, FID{findex2});
                    fname_savefile = fullfile(aas_getstudypath(aap), 'pairwise_maps', temp);
                    template_header.fname = fname_savefile;
                    template_header.dt = [spm_type('float32') 0];
                    spm_write_vol(template_header, RMAP);
                    individual_map_fnames{end+1} = fname_savefile; % for desc
                    
                end

            end % findex2
            
        end % findex1
            
        % save the average map(s)
 
        average_map_fname = fullfile(aas_getstudypath(aap),'ISC_average.nii');
        template_header.fname = average_map_fname;
        template_header.dt = [spm_type('float32') 0];
        RMAP = SUMMAP/file_count;

        % we always convert to Z before summing for averaging
        % if  convert_r_to_z just write
        % if ~convert_r_to_z, inverse transform and write
        
        if (convert_r_to_z == true)
            template_header.descrip = 'average ISC (z values)';
        else
            RMAP = tanh(RMAP);
            template_header.descrip = 'average ISC (r values)';
        end 

        spm_write_vol(template_header, RMAP);
        
        % desc the average map
        % optionally desc the individual maps
                        
        aap = aas_desc_outputs(aap,'average_intersubject_correlation_map', average_map_fname);
         
        if (save_individual_maps == true)
            aap = aas_desc_outputs(aap,'intersubject_correlation_maps', individual_map_fnames);
        end  
 
        % a summary correlation is saved but not desc'ed
        
        if (generate_summary_correlation_matrix == true)
            fname = fullfile(aas_getstudypath(aap),'summary_correlation_matrix.mat');
            save(fname,'summary_correlation_matrix_r','summary_correlation_matrix_p', 'summary_descriptors');
        end
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,true,sprintf('Unknown task %s',task));
        

end % switch

end









