function [aap,resp] = aamod_QA(aap, task)
%
% generate a QA summary of functional data
%
% Change History
%
% fall 2023 [MSJ] -- new
%

resp='';

switch task
	
    case 'report'
        
    case 'doit'		  
        
        nvoxels = aas_getsetting(aap,'nvoxels');
        carpetscale = aas_getsetting(aap,'carpetscale');

        outlier_filter = aas_getsetting(aap,'outlier_filter');
        outlier_threshold = aas_getsetting(aap,'outlier_threshold');
       
        imagedir = fullfile(aas_getstudypath(aap),'QA_images');
        if ~exist(imagedir,'dir');mkdir(imagedir);end
        
        datadir = fullfile(aas_getstudypath(aap),'QA_data');
        if ~exist(datadir,'dir');mkdir(datadir);end
        
        [ ~,subindex_list ] = aas_getN_bydomain(aap,'subject');

        summary_fname = fullfile(aas_getstudypath(aap),'data_summary.txt');

        if exist(summary_fname,'file');delete(summary_fname);end

        fid = fopen(summary_fname,'w'); % plaintext stats summary....
        
        if (fid < 0)
            aas_log(aap, true, 'Cannot open summary file for writing. Halting...');
        end
        
        group_outliers = []; % accumulate outliers for group outlier plot...

        for subject_index = subindex_list

            [ ~,all_sess ] = aas_getN_bydomain(aap, 'session', subject_index);
            session_list = intersect(all_sess, aap.acq_details.selected_sessions);           
            
            for session_index = session_list        

                fname = aas_getfiles_bystream(aap, subject_index, session_index, aas_getstreams(aap,'input',1)); 
        
                header = spm_vol(fname);
                data = spm_read_vols(header);
            
                [ r,c,s,f ] = size(data);
                data = reshape(data,r*c*s,f); % nvoxels x nframes
                                
                % stats

                nancount = sum(isnan(data(:)));
                zerocount = sum(data(:)==0);
                voxsize = diag(header(1).mat);
                voxsize = voxsize(1:3);

                % swap nan for zeros -- so masking can strip nan

                if (nancount > 0)
                    data(isnan(data)) = 0;
                end

                datamax = max(data(:));
                datamin = min(data(:));
                
                % extract voxel subset for processing / plotting
                
                selector = floor(linspace(1,r*c*s,nvoxels));
                data = data(selector,:);

                % we don't have a brainmask, so instead do a
                % quick and dirty implicit mask: we assume any
                % identically zero rows are outside of brain
                % this results in a reasonably good selection
                % especially if nvoxels is large-ish
                
                data = data(all(data~=0,2),:);  
  
                % removing DC offset improves joy and carpet plot color scaling...
                % (and has no effect on outlier identification)
                
                data = data - mean(data(:));
 
                % style points - open window in center of display
                % (this trick only works if running localsingle)
               
                if (strcmp(aap.options.wheretoprocess, 'localsingle'))
                    h = figure('Position',[0 0 1500 1000], 'Color', [1 1 1],'NumberTitle','off', 'Visible', 'off', 'MenuBar', 'none');
                    movegui(h, 'center');
                    set(h, 'Visible', 'on');
                else
                    h = figure('Position',[0 0 1500 1000],'Color', [1 1 1], 'NumberTitle','off','MenuBar','none');
                end

                subplot(3,3,[5 6 8 9]);

                % joy plot
                
                surf(data,log(abs(data)));
                axis tight
                xlabel('frame','FontSize',18);
                ylabel('voxel','FontSize',18);
                set(gca,'FontSize',18);
                shading interp;
                colormap gray;
                cm = colormap;
                colormap(cm.*cm.*cm); % suppresses midrange
                view(-25,60);

                % traditional carpet plot for comparison

                subplot(3,3,[4 7]);
                imagesc(data,[carpetscale(1)*max(data(:)) carpetscale(2)*max(data(:))]);
                xlabel('frame','FontSize',14);
                ylabel('voxel','FontSize',14);
                title('grayplot','FontSize',14);

                % average over voxels plot and outliers

                subplot(3,3,[1 2 3]);

                summary = mean(data,1,'omitnan');
                plot(summary,'k','LineWidth',2);
                a = axis;
                hold on;
                
                switch outlier_filter
                    
                    case 'anymedian'
                                        
                        outliers = isoutlier(data,'median',2, 'ThresholdFactor',outlier_threshold);
                        outliers = any(outliers,1);

                    case 'anymean'
                        
                        outliers = isoutlier(data,'mean',2, 'ThresholdFactor',outlier_threshold);
                        outliers = any(outliers,1);
                  
                    case 'mediansum'
                        
                        outliers = isoutlier(summary,'median', 'ThresholdFactor',outlier_threshold);                  
                        
                    case 'meansum'
                        
                        outliers = isoutlier(summary,'mean', 'ThresholdFactor',outlier_threshold);
 
                    case 'none'
                        
                        outliers = [];
                        
                end             
                
                outlier_index = find(outliers);
                if ~isempty(outlier_index)
                    plot(outlier_index,0.5*(a(3)+a(4)),'r+','LineWidth',2,'MarkerSize',12);
                end
                axis tight;
                
                % include stats in title
                
                ID = sprintf('%s-%s',aas_getsubjname(aap,subject_index),aas_getsessname(aap,session_index));

                titlestring = [strrep(ID,'_','-') '   dim: ' num2str([r c s f]) '   voxsize: ' num2str(voxsize')  '  min: ' num2str(datamin,3)  '  max: ' num2str(datamax,3)  '  nan: ' num2str(nancount) ' zeros: ' num2str(zerocount)];
                title(titlestring,'FontSize',18);
                
                if ~isempty(outlier_index)
                    odescriptor = sprintf('outliers (%s > %g)', outlier_filter, outlier_threshold);
                    legend({'global signal', odescriptor},'FontSize',16);
                else
                    legend('global signal','FontSize',16);
                end
                xlabel('frame','FontSize',18);
                ylabel('intensity','FontSize',18);

                fname_out = fullfile(imagedir,sprintf('%s.jpg', strrep(ID,'-','_')));

                set(h,'Renderer','opengl');
                set(findall(h,'Type','text'),'FontUnits','normalized');
                h.InvertHardcopy = 'off'; % otherwise matlab changes our nice black bg to white
                print(h, '-djpeg', '-r150', fname_out);
                
                close(h);
                
                % save the summary data

                fname_out = fullfile(datadir,sprintf('%s.mat', strrep(ID,'-','_')));
                save(fname_out,'data','summary');
                            
                % collect outliers for group outlier plot
                                
                if (isempty(group_outliers))
                    group_outliers(1,:) = outliers;
                else
                    % must protect against unequal number-of-frames
                    [~,ncols] = size(group_outliers);
                    if (ncols < length(outliers))
                        % this will zero extend ALL rows of group_outliers
                        % to new width so we can add this data to bottom
                        group_outliers(1,length(outliers)) = 0;
                    end
                    if (ncols > length(outliers))
                        % this will extend current outliers to fit
                        % group_outliers (we assume frame(1) are aligned)
                        outliers(ncols) = 0;
                    end
                    group_outliers(end+1,:) = outliers;
                end
                    
                % add stats to summary file
                
                fprintf(fid,'%s    outliers: %d\n', titlestring, sum(outliers));
              
            end % session
           
        end % subject
             
        fclose(fid); % summary_fname
        
        % group outlier plot
        % this is useful for detecting correlated artifacts across subjects
        
        % get a little better appearance if we size the window to fit data
               
        [nrows,ncols] = size(group_outliers);
        
        % don't need to make a group plot if there's just one subject
 
        if (nrows > 1)

            if (strcmp(aap.options.wheretoprocess, 'localsingle'))
                h2 = figure('Position',[0 0 5*ncols 50+20*nrows], 'Color', [1 1 1],'NumberTitle','off', 'Visible', 'off', 'MenuBar', 'none');
                movegui(h2, 'center');
                set(h2, 'Visible', 'on');
            else
                h2 = figure('Position',[0 0 5*ncols 50+20*nrows],'Color', [1 1 1], 'NumberTitle','off','MenuBar','none');
            end

            group_outliers(end+1,end+1) = 0; % pcolor needs a terminal dummy row/col
            pcolor(group_outliers);
            title('Outliers','FontSize',16);
            xlabel('frame','FontSize',14);
            ylabel('subject','FontSize',14);
            colormap gray;

            fname_out = fullfile(aas_getstudypath(aap),'group_outliers.jpg');
            set(h2,'Renderer','opengl');
            set(findall(h2,'Type','text'),'FontUnits','normalized');
            h2.InvertHardcopy = 'off';
            print(h2, '-djpeg', '-r150', fname_out);
            close(h2);
        
        end
        
        aap = aas_desc_outputs(aap, 'data_summary', summary_fname);

        % done!
               
    case 'checkrequirements'
        
    otherwise
        
        aas_log(aap,true,sprintf('Unknown task %s',task));
        
end % switch
end % function