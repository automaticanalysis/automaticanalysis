function [aap,resp] = aamod_DVARS(aap,task,subj,sess)
%
% generate epi QA metrics DVARS (= rms(diff(epi(:,:,:,t)))
% 

resp='';

switch task
	
    case 'domain'
		
        resp = 'session';
		
    case 'description'
		
        resp = 'compute wholebrain and, optionally, gm- and wm-specific DVARS';

    case 'report'
		    
    case 'doit'
		
		% your options are whole-brain DVARS (always generated), and
		% optionally grey and white matter DVARS (if masks are supplied)
					
		%------------------------------------------------------------------
		% WHOLE BRAIN 
		%------------------------------------------------------------------
		
		epi_file_names = aas_getfiles_bystream(aap, subj, sess, 'epi');
		epi_headers = spm_vol(epi_file_names);
		epi = spm_read_vols(epi_headers);
		
		[ row,col,slice,time ] = size(epi);
		
		epi = reshape(epi,row*col*slice,time);
				
		% quick-n-dirty implicit masking -- remove zero rows
		% (don't modify "epi" here! it may be needed later for gm
		% or wm extraction and removing rows will mess up the mask)
		
		temp = diff(epi,1,2);
		temp = temp(any(temp,2),:);
		wb = rms(temp);
				
		%------------------------------------------------------------------
		% GREY MATTER (if gm mask provided) 
		%------------------------------------------------------------------
	
		% sanity check -- the mask must match the epi. It won't if the
		% user didn't bother to reslice it in the tasklist (masks come
		% from segmentation and segmentation comes from the structural)
				
		if (aas_stream_has_contents(aap, subj, sess, 'native_greymask'))

			gm_mask_filename = aas_getfiles_bystream(aap, subj, 'native_greymask');
			mask_header = spm_vol(gm_mask_filename);
			representative_epi_header = epi_headers(1);

			if ( ~isequal(mask_header.dim, representative_epi_header.dim) || ...
					norm(mask_header.mat-representative_epi_header.mat)>0.01 )

				% fix or bail?

				if (aap.tasklist.currenttask.settings.resliceMaskIfNecessary)

					aas_log(aap, false, sprintf('\n%s: Reslicing gm mask to match epi...\n', mfilename));

					resliceOpts = [];
					resliceOpts.mask = false;
					resliceOpts.mean = false;
					resliceOpts.interp = 0;		% NN (because it's a MASK)
					resliceOpts.which = 1;		% don't reslice the first image
					resliceOpts.wrap = [0 0 0];
					resliceOpts.prefix = 'r';

					% need to pick off first frame of epi, otherwise
					% spm_reslice also reslices 2:end frames (which=1
					% applies to first IMAGE, not first file)
					
					spm_reslice({[representative_epi_header(1).fname ',1'] gm_mask_filename}, resliceOpts);
				
					% sneaklily replace gm_mask_filename so we spm load resliced version

					[ p,n,e ] = fileparts(gm_mask_filename);
					gm_mask_filename = fullfile(p,[resliceOpts.prefix n e]);

				else
					aas_log(aap, true, sprintf('\nGM Mask incompatible with epi. Reslice in tasklist prior to using %s\n', mfilename));
				end

			end

			gm_mask = spm_read_vols(spm_vol(gm_mask_filename));
			
			% ~~ is a clever trick to make sure mask is 0/1 (not 0/something)
			
			gm = rms( diff(epi(~~gm_mask(:),:),1,2) );	% grey matter DVARS
						
		else
			gm = [];
		end
		
		%------------------------------------------------------------------
		% WHITE MATTER (if wm mask provided)
		%------------------------------------------------------------------
				
		if (aas_stream_has_contents(aap, subj, sess, 'native_whitemask'))
			
			wm_mask_filename = aas_getfiles_bystream(aap, subj, 'native_whitemask');
			mask_header = spm_vol(wm_mask_filename);
			representative_epi_header = epi_headers(1);
	
			if ( ~isequal(mask_header.dim, representative_epi_header.dim) || ...
					norm(mask_header.mat-representative_epi_header.mat)>0.01 )

				% fix or bail?

				if (aap.tasklist.currenttask.settings.resliceMaskIfNecessary)

					aas_log(aap, false, sprintf('\n%s: Reslicing wm mask to match epi...\n', mfilename));

					resliceOpts = [];
					resliceOpts.mask = false;
					resliceOpts.mean = false;
					resliceOpts.interp = 0;		% NN (because it's a MASK)
					resliceOpts.which = 1;		% don't reslice the first image
					resliceOpts.wrap = [0 0 0];
					resliceOpts.prefix = 'r';

					spm_reslice({[representative_epi_header(1).fname ',1'] wm_mask_filename}, resliceOpts);

					[ p,n,e ] = fileparts(wm_mask_filename);
					wm_mask_filename = fullfile(p,[resliceOpts.prefix n e]);

				else
					aas_log(aap, true, sprintf('\nWM Mask incompatible with epi. Reslice in tasklist prior to using %s\n', mfilename));
				end

			end

			wm_mask = spm_read_vols(spm_vol(wm_mask_filename));
			
			wm = rms( diff(epi(~~wm_mask(:),:),1,2) );	% white matter DVARS
			
		else
			wm = [];
		end
									
		% desc the outputs
		%
		% this will create DVARS.wb, DVARS.gm, and DVARS.wm
		% on subsequent DVARS = load('DVARS.mat')
		
		outfile = fullfile(aas_getsesspath(aap,subj,sess), 'DVARS.mat'); 
		save(outfile,'wb','gm','wm');
		aap = aas_desc_outputs(aap, subj, sess, 'DVARS', outfile);

		
							
end