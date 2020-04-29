function [aap,resp] = aamod_greyplot(aap,task,subj,sess)
%
% generate epi greyplot and also four user-supplied metrics
%
% plots are saved as jpegs. additionally, the data from the four
% metric plots are saved as the outputstream metric_data. The struct is
% organized metric_data.fieldname (one for each defined metric), where
% fieldname is taken from the "fieldname" field in the xml header entry for
% the metric
%
% The greyplot will be restricted to greymatter voxels if a gm mask is
% provided, whitematter voxels if a wm mask is provided and a 50/50 split if
% gm and wm voxel if both are provided. If no masks are provided, the
% greyplot is based on whole-brain voxels. The number of voxels include in
% the plot is set by numberOfVoxelsToPlot in the header -- this is used as
% the upper limit of a linear selection of all voxels (whole-brain,
% wm-only, gm-only, or wm+gm, depending on the masks provided)
%
%
% Reference:
%
% 	Jonathan D. Power
% 	A simple but useful way to assess fMRI scan qualities
% 	NeuroImage 154, 150-158, (2017) 
% 	https://doi.org/10.1016/j.neuroimage.2016.08.009
%

resp='';

switch task
	
    case 'domain'
		
        resp = 'session';
		
    case 'description'
		
        resp = 'generate and save epi greyplot and related metrics';

    case 'report'
		
		sesspath = aas_getsesspath(aap, subj, sess);
		greyplot = fullfile(sesspath,'greyplot.jpg');
		aap = aas_report_add(aap, subj, '<table><tr><td>');
		aap = aas_report_addimage(aap, subj, greyplot);
		aap = aas_report_add(aap, subj, '</td></tr></table>');
    
    case 'doit'

		% required stream input 
				
		epi_file_names = aas_getfiles_bystream(aap, subj, sess, 'epi');
		epi_headers = spm_vol(epi_file_names);
		epi = spm_read_vols(epi_headers);
		
		[ nrow,ncol,nslice,nframes ] = size(epi);
		voxel_count = nrow*ncol*nslice;

		% smoothing epi gives a better looking greyplot but we need to do this before
		% we destroy the 3D structure. Ergo, filter each frame, then unroll into
		% a 2D (voxels x time) array wb_epi ("wholebrain epi")
		
		% update -- use spm_smooth so we don't require image processing toolbox

		wb_epi = zeros(voxel_count,nframes);
		temp = zeros(size(epi(:,:,:,1)));
		for frame = 1:nframes
			spm_smooth(epi(:,:,:,frame) ,temp, 3);
			wb_epi(:,frame) = temp(:);
		end
		wb_epi = detrend(wb_epi')';
			
		temp = aas_getfiles_bystream(aap, subj, sess, 'GLOBALMEAN');
		temp = load(temp);
		GLOBALMEAN = temp.GLOBALMEAN;
        
		temp = aas_getfiles_bystream(aap, subj, sess, 'FD');
		temp = load(temp);
		FD = temp.FD;

		temp = aas_getfiles_bystream(aap, subj, sess, 'DVARS');
		DVARS = load(temp); % creates DVARS.wb, DVARS,wm, DVARS.gm

		% ----------------------------------------------------------------
		% extract gm voxels (optional)
		% ----------------------------------------------------------------
		
		if (aas_stream_has_contents(aap, 'native_greymask'))
			
			% sanity check -- the mask must match the epi. It won't if the
			% user didn't bother to reslice it in the tasklist (masks come
			% from segmentation and segmentation comes from the structural)

			gm_mask_filename = aas_getfiles_bystream(aap, subj, 'native_greymask');
			mask_header = spm_vol(gm_mask_filename);
			
			% we test/reslice against mean epi
			
			meanepi_filename = aas_getfiles_bystream(aap, subj, 'meanepi');
			meanepi_header = spm_vol(meanepi_filename);

			if ( ~isequal(mask_header.dim, meanepi_header.dim) || ...
					norm(mask_header.mat-meanepi_header.mat)>0.01 )

				% fix or bail?

				if (aap.tasklist.currenttask.settings.resliceMaskIfNecessary)

					aas_log(aap, false, sprintf('\n%s: Reslicing gm mask to match wb_epi...\n', mfilename));

					resliceOpts = [];
					resliceOpts.mask = false;
					resliceOpts.mean = false;
					resliceOpts.interp = 0;		% NN (because it's a MASK)
					resliceOpts.which = 1;		% DON'T reslice the first image
					resliceOpts.wrap = [0 0 0];	% this is everywhere in aa even though it should be [1 1 0] for MRI
					resliceOpts.prefix = 'r';
					
					% meanepi_header.fname is 4D -- spm-reslice
					% is reslicing vols 2-n in addition to structural!				
					
					spm_reslice({meanepi_header.fname gm_mask_filename}, resliceOpts);

					% sneaklily replace gm_mask_filename so we spm load resliced version
					
					[ p,n,e ] = fileparts(gm_mask_filename);
					gm_mask_filename = fullfile(p,[resliceOpts.prefix n e]);

				else
					aas_log(aap, true, sprintf('\nWM Mask incompatible with epi. Reslice in tasklist prior to using %s\n', mfilename));
				end

			end
	
			mask_header = spm_vol(gm_mask_filename);
			gm_mask = spm_read_vols(mask_header);
	
			% sanity check -- epi needs to match the mask (otherwise the voxel extraction will crash because
			% of an array size mismatch). This can happen if you run this module after normalization (because 
			% aa will feed it the most recent epi, which is the normed epi). Either explicitly pick an unnormed
			% epi stream (see aa documentation) or put this module *before* normalization in the tasklist
			
			if ( ~isequal(mask_header.dim, epi_headers(1).dim) || ...
					norm(mask_header.mat-epi_headers(1).mat)>0.01 )
						
				aas_log(aap, true, sprintf('\n%s: epi is not in native space. Run greyplot BEFORE normalization...\n', mfilename));
			
			end

			% ~~ is a clever trick to make sure mask is 0/1 (not 0/something)

			gm_masked_epi = wb_epi(~~gm_mask(:),:);			
				
		else
			gm_masked_epi = [];
		end
		
		[ gm_voxel_count,~ ] = size(gm_masked_epi);
		
		% ----------------------------------------------------------------
		% extract wm voxels (optional)
		% ----------------------------------------------------------------
		
		if (aas_stream_has_contents(aap, 'native_whitemask'))
			
			wm_mask_filename = aas_getfiles_bystream(aap, subj, 'native_whitemask');
			mask_header = spm_vol(wm_mask_filename);
			meanepi_filename = aas_getfiles_bystream(aap, subj, 'meanepi');
			meanepi_header = spm_vol(meanepi_filename);

			if ( ~isequal(mask_header.dim, meanepi_header.dim) || ...
					norm(mask_header.mat-meanepi_header.mat)>0.01 )

				if (aap.tasklist.currenttask.settings.resliceMaskIfNecessary)

					aas_log(aap, false, sprintf('\n%s: Reslicing gm mask to match wb_epi...\n', mfilename));

					resliceOpts = [];
					resliceOpts.mask = false;
					resliceOpts.mean = false;
					resliceOpts.interp = 0;
					resliceOpts.which = 1;
					resliceOpts.wrap = [0 0 0];
					resliceOpts.prefix = 'r';

					spm_reslice({meanepi_header.fname wm_mask_filename}, resliceOpts);

					[ p,n,e ] = fileparts(wm_mask_filename);
					wm_mask_filename = fullfile(p,[resliceOpts.prefix n e]);

				else
					aas_log(aap, true, sprintf('\nWM Mask incompatible with wb_epi. Reslice in tasklist prior to using %s\n', mfilename));
				end

			end
					
			mask_header = spm_vol(wm_mask_filename);
			wm_mask = spm_read_vols(mask_header);
				
			if ( ~isequal(mask_header.dim, epi_headers(1).dim) || ...
					norm(mask_header.mat-epi_headers(1).mat)>0.01 )
						
				aas_log(aap, true, sprintf('\n%s: epi is not in native space. Run greyplot BEFORE normalization...\n', mfilename));
			
			end

			% ~~ is a clever trick to make sure mask is 0/1 (not 0/something)

			wm_masked_epi = wb_epi(~~wm_mask(:),:);	
		
		else
			wm_masked_epi = [];
		end
		
		[ wm_voxel_count,~ ] = size(wm_masked_epi);		
			
		% ----------------------------------------------------------------
		% plot greyscale-related QA metrics
		% ----------------------------------------------------------------
		
		% the movegui trick here centers the window, which looks cool
		% but it breaks plotting if running headless (i.e. cluster)
		
		if (strcmp(aap.options.wheretoprocess, 'localsingle'))
			figh = figure('Color',[1 1 1],'Position',[0 0 750 960],'Visible', 'off');
			movegui(figh,'center');
			set(figh,'Visible', 'on');
		else
			figh = figure('Color',[1 1 1],'Position',[0 0 750 960]);
		end
		
		metric_plots = aap.tasklist.currenttask.settings.metric_plots;
		plot_count = numel(metric_plots.plot);
				
		for index = 1:plot_count
		
			% the subplots go three across (for spacing)
			% and the greyplot gets 2
			
			m=3*(index-1)+1;n=3*(index-1)+2;
			
			subplot(plot_count+2,3,[m,n]);

			this_plot_data = eval(metric_plots.plot(index).eval_string);
			this_plot_data = this_plot_data(:); % force col vec
			
			plot(this_plot_data);
			axis tight;

			if (metric_plots.plot(index).label3sig)
				threesig = nanmean(this_plot_data) + 3*mad(this_plot_data);
				h = refline(0,threesig);
				h.Color = 'r';
			end

			ylabel(metric_plots.plot(index).yaxis_label);

			% add a histogram
			ax=gca;p=ax.Position;
			pvector=[ p(1)+p(3) p(2) 0.2 p(4)];
			[ counts,centers ] = hist(this_plot_data,50);
			subplot('position', pvector);
			barh(centers, counts, 'FaceColor',[0.6 0.6 0.6],'EdgeColor',[0.6 0.6 0.6]);
			axis( [0 Inf -Inf Inf] );
			axis off;
			
			% while we're here, save the data for later stream output
			
			% we save the eval string so subsequent modules can verify how
			% the data was generated. We need to save the data in a struct
			% in case they are different lengths. However, eval_string may 
			% not be a valid struct field name because of whitespace
			% and operators. Ergo, we require a fieldname for the data:
			
			metric_data.(metric_plots.plot(index).fieldname) = this_plot_data;
			
		end
		
		subplot(plot_count+2,3,[1,2]);
		title(strrep([ aas_getsubjname(aap,subj) '/' aas_getsessname(aap,sess) ],'_','-'));

		% ----------------------------------------------------------------
		% add the greyplot
		% ----------------------------------------------------------------
		
		% if no mask provided, use whole-brain epi
		% if one mask provided, use masked
		% if both masks provided, use half of both
		
		numberOfVoxelsToPlot = aap.tasklist.currenttask.settings.numberOfVoxelsToPlot;
		
		if (~isempty(gm_masked_epi) && ~isempty(wm_masked_epi))
			numberOfVoxelsToPlot = round(numberOfVoxelsToPlot/2);
		end
		
		smooth_img = [];
		
		% conv2(ones(1,2)... is just a cheap 2D smoothing filter to clean up the map a bit
		
		if (isempty(gm_masked_epi) && isempty(wm_masked_epi))
			voxel_subset = round(linspace(1, voxel_count, numberOfVoxelsToPlot));
			smooth_img = conv2(ones(1,2),wb_epi(voxel_subset,:));
		end
			
		if (~isempty(gm_masked_epi))
			voxel_subset = round(linspace(1, gm_voxel_count, numberOfVoxelsToPlot));
			smooth_img = [ smooth_img ;  conv2(ones(1,2),gm_masked_epi(voxel_subset,:))];
		end
		
		if (~isempty(wm_masked_epi))
			voxel_subset = round(linspace(1, wm_voxel_count, numberOfVoxelsToPlot));
			smooth_img = [ smooth_img ;  conv2(ones(1,2),wm_masked_epi(voxel_subset,:))];
		end
		
		% the subplots go three across (for spacing)
		% and the greyplot occupies the last two rows
		
		subplot(plot_count+2,3,[3*plot_count+1 3*plot_count+2 3*(plot_count+1)+1 3*(plot_count+1)+2]);
		
		imagesc(smooth_img);
		xlabel('frame');ylabel('voxel');
		colormap(gray); % curse you Mathworks for your rebellious American spelling!
		caxis(0.8*caxis); % decomment to pump up the contrast a bit...
		ax=gca;p=ax.Position;
		pvectors=[ p(1)+p(3)+0.02 p(2) 0.02 p(4) ];
		hb=colorbar;
		set(hb, 'Units', 'normalized', 'position', pvectors);

		% save to jpeg
		
		fname = fullfile(aas_getsesspath(aap,subj,sess), 'greyplot.jpg'); 
		set(figh,'Renderer','opengl');
		set(findall(figh,'Type','text'),'FontUnits','normalized');
		print(figh, '-djpeg', '-r150', fname);
		close(figh);

			
		% ----------------------------------------------------------------
		% desc the metric data
		% ----------------------------------------------------------------
		
		fname = fullfile(aas_getsesspath(aap,subj,sess), 'metric_data.mat'); 
		save(fname,'metric_data');
		aap = aas_desc_outputs(aap, subj, sess, 'metric_data', fname);


							
end