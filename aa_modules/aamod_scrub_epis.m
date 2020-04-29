function [aap,resp] = aamod_scrub_epis(aap, task, subj, sess)
%
% aamod_scrub_epis
%
%	scrub epi volumes based on criteria supplied
%
%	NB: this module generates a list of delta regressors to "scrub" data
%	from the GLM -- it does not literally delete the frames from the nii
%
% Revision History
%
% winter 2018 [MSJ] -- absorb aamod_listspike 
% spring 2018 [MSJ] -- new
%

resp='';

switch task
	
    case 'report'
        
    case 'doit'
	
		% initialization stuff
		
		explicit_scrublist = [ ];
		motion_scrublist = [ ];
		master_scrublist = [ ];
		
		% stream input ----------------------------------------------------------
		
		session_path = aas_getsesspath(aap, subj, sess);
		
		% get the number of volumes (frames) in the epi
		
		epi_file_list = aas_getimages_bystream(aap, subj, sess, 'epi');
		numvol = size(epi_file_list, 1);
		
		% if 4D epi (i.e. one file) we need to extract numvol from header
		
		if (numvol == 1)
			numvol = numel(spm_vol(epi_file_list));
		end
		
		% generate list of frames to scrub  -------------------------------------	
				
		% 1) begin with frames from explicit_framelist if it's defined
	
		if (~isempty(aap.tasklist.currenttask.settings.explicit_framelist))
			explicit_framelist = aap.tasklist.currenttask.settings.explicit_framelist;
			explicit_scrublist = eval(explicit_framelist);
		end
		

		% 2) add a metric-based scrublist using scrub_criteria if the var is not empty

		if (~isempty(aap.tasklist.currenttask.settings.scrub_criteria))
			
			% load whatever input streams are defined

			% is there a better way to do this? we'd like to only load the
			% streams the scrub_criteria is going to use, but that seems
			% prone to error...
			
			if (aas_stream_has_contents(aap, subj, sess, 'realignment_parameter'))
				rp = aas_getfiles_bystream(aap, subj, sess,'realignment_parameter');
				rp = spm_load(rp);
				rpdiff = [zeros(1,6); diff(rp)];
				dXYZ = rpdiff(:,1:3);
				dRPJ = rpdiff(:,4:6);
			end
 
			if (aas_stream_has_contents(aap, subj, sess, 'FD'))
				temp = aas_getfiles_bystream(aap, subj, sess, 'FD');
				temp = load(temp);
				FD = temp.FD;
			end

			if (aas_stream_has_contents(aap, subj, sess, 'DVARS'))
				temp = aas_getfiles_bystream(aap, subj, sess, 'DVARS');
				DVARS = load(temp); % creats DVARS.wb, DVARS.wm, DVARS.gm
			end

			if (aas_stream_has_contents(aap, subj, sess, 'tsdiffana'))
				temp = aas_getfiles_bystream(aap, subj, sess, 'tsdiffana');
				tsdiffana = load(temp);
				SV = tsdiffana.qa.global.diff / mean(tsdiffana.qa.global.mean);
				GLOBALS = tsdiffana.qa.global.mean;
 				% this was definedin listspikes - CHECK
				TM = tsdiffana.qa.global.diff / mean(tsdiffana.qa.global.mean).^2;
			end

			if (aas_stream_has_contents(aap, subj, sess, 'metric_data'))
				temp = aas_getfiles_bystream(aap, subj, sess, 'metric_data');
				temp = load(temp);
				metric_data = temp.metric_data;
			end
			
			if (aas_stream_has_contents(aap, 'metric_thresholds'))
				temp = aas_getfiles_bystream(aap, 'metric_thresholds');
				temp = load(temp);
				metric_thresholds = temp.metric_thresholds;
			end
			
			% unfortunately, there's lots of opportunity for error here and not much we can do about it...
			
			try
				
				temp = eval(aap.tasklist.currenttask.settings.scrub_criteria);
				motion_scrublist = find(temp);
				
			catch
				
				aas_log(aap, true, sprintf('%s: eval of scrub criterion %s failed.\n', mfilename, aap.tasklist.currenttask.settings.scrub_criteria));

			end
			
		end
						
		%  "union" will sort lists and remove duplicates
		
		master_scrublist = union(explicit_scrublist, motion_scrublist);
		
		% we may extend scrubbing a fixed number of frames before and
		% after any identified frames targeted for scrubbing - the
		% range is specified by the prekillbox and postkillbox 
		% parameters (0 = don't extend)
		
		% easiest if we convert to indicator variables to do this:
		
		master_scrub_indicator = zeros(numvol,1);
		master_scrub_indicator(master_scrublist) = 1;

		if (~isempty(aap.tasklist.currenttask.settings.prekillbox))
			prekillbox = aap.tasklist.currenttask.settings.prekillbox;
		else
			prekillbox = 0;
		end
		
		if (prekillbox > 0)
			for index=1:length(master_scrub_indicator)
				if (master_scrub_indicator(index))
					istart = index - prekillbox;
					istart = max(istart,1);
					master_scrub_indicator(istart:index) = 1;
				end
			end
		end
				
		if (~isempty(aap.tasklist.currenttask.settings.postkillbox))
			postkillbox = aap.tasklist.currenttask.settings.postkillbox;
		else
			postkillbox = 0;
		end
		
		if (postkillbox > 0)
			for index=length(master_scrub_indicator):-1:1
				if (master_scrub_indicator(index))
					iend = index + postkillbox;
					iend = min(iend,length(master_scrub_indicator));
					master_scrub_indicator(index:iend) = 1;
				end
			end
		end		
		
		% process indicator list into scrub and keeper lists
							
		scrublist = []; keeperlist = [];
		for index = 1:length(master_scrub_indicator)
			if (master_scrub_indicator(index))
				scrublist = [ scrublist index ];
			else
				keeperlist = [ keeperlist index ];
			end
		end
		
  
        %  diagnostic images  ---------------------------------------------
				
		if (strcmp(aap.options.wheretoprocess, 'localsingle'))
			% centering trick only works if not cluster
			h = figure('Position',[0 0 500 100], 'Visible', 'off', 'MenuBar', 'none');
			movegui(h, 'center');
			set(h, 'Visible', 'on');
		else
			h = figure('Position',[0 0 500 100],'MenuBar','none');
		end
		master_scrub_indicator = [ master_scrub_indicator ; 0 ]; % for pcolor weirdness
		master_scrub_indicator = [ master_scrub_indicator' ; master_scrub_indicator' ];
		pcolor(master_scrub_indicator);
		axis off; colormap('flag');
		title(strrep([ aas_getsubjname(aap,subj) '/' aas_getsessname(aap,sess) ],'_','-'));
		set(h,'Renderer','opengl');
		set(findall(h,'Type','text'),'FontUnits','normalized');
		fname = fullfile(session_path, 'keeplist.jpg');
		print(h, '-djpeg', '-r150', fname);
		close(h);

		% (there's more diag images from aamod_listspikes we could include here ) 
		
        %  desc -----------------------------------------------------------
		
		% save the scrub and keeper lists. These are just a list of integers.
		% save( ) converts these to double, so we use dlmwrite( ) here which
		% doesn't. Either works with load( ) if you need the data later...
		      
        scrublist_fname = fullfile(session_path, 'scrublist.txt');
		dlmwrite(scrublist_fname, scrublist);
		aap = aas_desc_outputs(aap, subj, sess, 'scrublist', scrublist_fname);
				
		keeperlist_fname = fullfile(session_path, 'keeperlist.txt');
		dlmwrite(keeperlist_fname, keeperlist);
		aap = aas_desc_outputs(aap, subj, sess, 'keeperlist', keeperlist_fname);
        
        % aamod_firstlevel_model uses the stream name "listspikes" (not
        % "scrublist") so we have to save a copy of the scrub info to this 
		% streamname if we expect to use it in the GLM
		
		% also, have to kludge our scrublist into TSspikes and Mspikes --
		% aas_firstlevel_model_nuisance expects these to be column vectors
		% and does a union on them, so here we set both equal to transpose
		% of the scrublist (which gets generated above as a row vec)
		
		% also be advised aamod_firstlevel_model won't use this unless
		% you truthy the includespikes option
					
		TSspikes=scrublist'; Mspikes=scrublist';
		listspikes_fname = fullfile(session_path, 'listspikes.mat');
		save(listspikes_fname, 'TSspikes', 'Mspikes');
        aap = aas_desc_outputs(aap, subj, sess, 'listspikes', listspikes_fname);

		
    case 'checkrequirements'
        
    otherwise
        aas_log(aap, 1, sprintf('%s: Unknown task %s', mfilename, task));
		
		
end