function [aap,resp] = aamod_GLOBALMEAN(aap, task, subj, sess)
%
% extract globalmean signal

resp='';

switch task
	
    case 'report'

	case 'doit'
		
		sesspath = aas_getsesspath(aap,subj,sess);
		
 		epi_file_list = aas_getimages_bystream(aap, subj, sess, 'epi');	
		numvol = size(epi_file_list, 1);
		epi_header = spm_vol(epi_file_list);
			
		template_mask_fname = aap.tasklist.currenttask.settings.template_mask;
		
		if ~isempty(template_mask_fname)
			
			if ~exist(template_mask_fname, 'file')
				aas_log(aap, true, sprintf('\n%s: Cannot find template mask %s Exiting...\n', mfilename, template_mask_fname));
			end
			
			% make a local copy of the template -- this is a little wasteful
			% but otherwise the rWawa gets created where the template lives
			
			[ ~,n,e ] = fileparts(template_mask_fname);
			tempname = fullfile(sesspath,[n e]);
			copyfile(template_mask_fname,tempname);
			template_mask_fname = tempname;
		
			mask_header = spm_vol(template_mask_fname);

			if (~isequal(mask_header.dim, epi_header(1).dim) ||	norm(mask_header.mat-epi_header(1).mat)>0.01)

				aas_log(aap, false, sprintf('\n%s: Reslicing mask to match epi...\n', mfilename));

				resliceOpts = [];
				resliceOpts.mask = false;
				resliceOpts.mean = false;
				resliceOpts.interp = 0;		% NN (because it's a MASK)
				resliceOpts.which = 1;		% DON'T reslice the first image
				resliceOpts.wrap = [0 0 0];	% this is everywhere in aa even though it should be [1 1 0] for MRI
				resliceOpts.prefix = 'r';

				% note selection of frame #1 (otherwise spm_reslices frames
				% 2-end in addition to reslicing templace_mask_fname
				
				spm_reslice({[epi_header(1).fname ',1'] template_mask_fname}, resliceOpts);

				% sneaklily overwrite template_mask_fname so we spm load resliced version
				
				[ p,n,e ] = fileparts(template_mask_fname);
				template_mask_fname = fullfile(p,[resliceOpts.prefix n e]);
				
				mask_header = spm_vol(template_mask_fname);
				
				% save a diagnostic image
				save_three_ortho_jpgs(template_mask_fname, fullfile(sesspath,'diagnostic_brainmask'));

			end
		
			mask = spm_read_vols(mask_header);
			have_explicit_mask = 1;
			
			% also turn off implicit mask, regardless of setting
		
			implicit_mask	= '';
			
		else			   
		
			% optional renameable explicit mask is first entry in inputstream
			% -- if defined, check the size and orientation

			first_inputstream_struct = aap.tasklist.currenttask.inputstreams(1).stream{1};

			if (aas_stream_has_contents(aap, first_inputstream_struct.CONTENT))
				try
					% mask_fname = aas_getfiles_bystream(aap, subj, sess, first_inputstream_struct.CONTENT);
					% the mask is created at subject level, not session level
					mask_fname = aas_getfiles_bystream(aap, subj, first_inputstream_struct.CONTENT);
					mask_header = spm_vol(mask_fname);
					% sanity check the mask against the data
					if (~isequal(mask_header.dim, epi_header(1).dim))
						aas_log(aap, true, sprintf('\n%s: Mask and data are different size. Exiting...\n', mfilename));
					end
					if (norm(mask_header.mat-epi_header(1).mat)>0.01)
						aas_log(aap, true, sprintf('\n%s: Mask and data have different orientations. Exiting...\n', mfilename));
					end
					mask = spm_read_vols(spm_vol(mask_header));
					have_explicit_mask = 1;
				catch
					aas_log(aap, true, sprintf('\n%s: Error reading mask %s. Exiting...\n', mfilename, first_inputstream_struct.CONTENT));
				end
			else
				have_explicit_mask = 0;
			end

			% implict_mask is just a number. Leave empty for "no implicit masking"

			implicit_mask = aap.tasklist.currenttask.settings.implicit_mask;
		
		end

		% process
		
		if (numvol == 1)
		
			% single 4D epi file
			
			V = spm_read_vols(epi_header);
			GLOBALMEAN = zeros(size(V,4),1);
			
			for t = 1:size(V,4)
				temp = V(:,:,:,t);
				if (have_explicit_mask) temp(~mask) = []; end
				temp = temp(:);
				if ~isempty(implicit_mask)
					temp(temp<implicit_mask) = [];
				end
				GLOBALMEAN(t) = mean(temp);
			end
			
		else
			
			% multiple 3D epi files
			
			GLOBALMEAN = zeros(numvol,1);
			
			for t = 1:numvol
		        temp = spm_read_vols(epi_header(t));
				if (have_explicit_mask) temp(~mask) = []; end
				temp = temp(:);
				if ~isempty(implicit_mask)
					temp(temp<implicit_mask) = [];
				end
				GLOBALMEAN(t) = mean(temp);
			end
			
		end	
		
		% describe output

		fname = fullfile(sesspath, 'GLOBALMEAN.mat');
		save(fname,'GLOBALMEAN');
		aap=aas_desc_outputs(aap, subj, sess, 'GLOBALMEAN', fname);
		
		% lets plot a diagnostic
				
		h = figure;
		plot(GLOBALMEAN);
		title('global signal');
		grid on; axis tight;
		set(h,'Renderer','opengl');
		% workaround for font rescaling weirdness
		set(findall(h,'Type','text'),'FontUnits','normalized');
		fname = fullfile(sesspath, 'diagnostic_globalsignal.jpg');
		print(h, '-djpeg', '-r150', fname);
		close(h);	

		% done!
 
	case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end

end
