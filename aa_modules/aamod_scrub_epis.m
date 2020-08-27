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
% 08/2020 [MSJ] - fix possible PCT race error
% summer 2020 [MSJ] -- added task timing plot
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
                
                % there's a weird random bug when using the PCT that a_g_b
                % sometimes returns an empty file for metric_thresholds
                % even though it's copied to the working directory properly.
                % Maybe a race condition? Anyway, pause and reload seems to
                % workaround the problem:
                
                if (~exist(temp,'file'))

                    aas_log(aap, false, sprintf('***** CANNOT FIND %s. Hail Mary pause ****** \n', temp));
                    unix('touch /Users/peellelab/SCRUB_EPI_FILE_FAIL.txt');
                    pause(1);
                    temp = aas_getfiles_bystream(aap, 'metric_thresholds');

                end        
                
				temp = load(temp);
				metric_thresholds = temp.metric_thresholds;
                
            end
  
            if (aas_stream_has_contents(aap, 'GLOBALMEAN'))
                temp = aas_getfiles_bystream(aap, subj, sess, 'GLOBALMEAN');
                temp = load(temp);
                GLOBALMEAN = temp.GLOBALMEAN;
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
  
        %  diagnostic bargraphs  -------------------------------------------------------------------------------------------
				        
        if (strcmp(aap.options.wheretoprocess, 'localsingle'))
            % centering trick only works if not cluster
            hf = figure('Position',[0 0 900 200], 'Color', [1 1 1],'NumberTitle','off', 'Visible', 'off', 'MenuBar', 'none');
            movegui(hf, 'center');
            set(hf, 'Visible', 'on');
		else
			hf = figure('Position',[0 0 900 200],'Color', [1 1 1], 'NumberTitle','off','MenuBar','none');
        end
        
        % bargraph settings
        
        bw = 0.95; % barwidth (0.8 = default; 1 = no overlap)
        bg = ones(size(master_scrub_indicator));
        
        % get TR from the epi header to convert event timing to frame #
        %
        % note we also need to know if timing is in secs or frames and
        % that info is only defined in SPM.xBF.UNITS which we can't access. 
        % Ergo, we added a setting to the header.
        
        if strcmp(aap.tasklist.currenttask.settings.xBFUNITS,'secs')
            DICOMHEADERS = load(aas_getfiles_bystream(aap,subj,sess,'epi_dicom_header'));
            TR = DICOMHEADERS.DICOMHEADERS{1}.volumeTR;
        else
            TR = 1;
        end
        
        % if there is no model defined, t_t_d will return empty event_names

        [ event_names, frameseries ] = task_timing_data(aap, subj, sess, numvol, TR);


        if (isempty(event_names))

            bh = bar(bg,bw);
            set(bh,'FaceColor',[0.7 0 0]);
            hold on;
            bar(master_scrub_indicator,bw,'k');
            axis tight; axis off;
            title(strrep([ aas_getsubjname(aap,subj) '/' aas_getsessname(aap,sess) ],'_','-'));

        else

            subplot(2,1,1);
            bh = bar(bg,bw);
            set(bh,'FaceColor',[0.7 0 0]);
            hold on;
            bar(master_scrub_indicator,bw,'k');
            axis tight; axis off;
            title(strrep([ aas_getsubjname(aap,subj) '/' aas_getsessname(aap,sess) ],'_','-'),'FontName','Helvetica', 'FontSize', 12);

            subplot(2,1,2)

            clist = [ 'r', 'g', 'b', 'c', 'm', 'y', 'k' ];
            cvals = 0.8 * [ 1 0 0 ; 0 1 0; 0 0 1; 0 1 1 ; 1 0 1; 1 1 0; 1 1 1 ];
            title_string = cell(numel(event_names),1);

            for index = 1:numel(event_names)
                bh = bar(frameseries(index,:),bw,clist(mod(index,length(clist))));
                set(bh,'FaceColor',cvals(index,:));
                taskloss = round(100 * dot(frameseries(index,:),master_scrub_indicator) / sum(frameseries(index,:)));
                title_string{index} = sprintf('%s(%s) %d%%     ', strrep(event_names{index},'_','-'), clist(mod(index,length(clist))), taskloss);
                hold on;
            end

            set(bh(1),'ShowBaseLine','off');
            axis tight;
            % this works better than axis off
            bh.Parent.YTick = []; bh.Parent.XTick = [];

            
            fsize = 12;
            if (numel(event_names) > 3); fsize=10; end
            if (numel(event_names) > 5); fsize=8; end
                
            title(sprintf('%s',title_string{:}),'FontName','Helvetica', 'FontSize', fsize);

        end    
        
   		set(hf,'Renderer','opengl');
		set(findall(hf,'Type','text'),'FontUnits','normalized');
		fname = fullfile(session_path, 'keeplist.jpg');
		print(hf, '-djpeg', '-r150', fname);
		close(hf);      
        		
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

end



% ---------------------------------------------------------------------------------------
% task_timing_data
% ---------------------------------------------------------------------------------------

function [ event_names, frameseries ] = task_timing_data(aap, subj, sess, nframes, TR)

% generate a task timing plot (TTP) for subject subj / session sess

% extract task timing data from aap struct for subject subj / session sess

% a "frameseries" is just a time series expressed as binary in frames
% (e.g., [ 4 5 6 ] => [ 0 0 0 1 1 1 0 0 0 0], assuming 10 frames)

event_names = [];
frameseries = [];

% must have firstlevel model settings -- bail if no

if (~isfield(aap.tasksettings,'aamod_firstlevel_model')) return; end
if (~isfield(aap.tasksettings.aamod_firstlevel_model, 'model')) return; end
    
model = aap.tasksettings.aamod_firstlevel_model.model;   

% we can't assume model struct is in the same order as subjects and
% sessions (i.e. model(subj) may not be correct) -- it depends on
% how aas_addevents were ordered -- so loop over all entries 
% and process matches (guess: there will only be one match because
% addevent collects all info for one subj/sess in one entry)

subject_name = aap.acq_details.subjects(subj).subjname;
session_name = aas_getsessname(aap,sess);

for mindex = 1:numel(model)
        
    if (strcmp(model(mindex).subject,subject_name) && strcmp(model(mindex).session, session_name))
        
        % we have a match
   
        event_list = model(mindex).event;

        % event_list will contain all event types (e.g., finger, foot, lips)
        % for this subject/session -- generate a frameseries for each one
        
        nevents = numel(event_list);
        
        event_names = cell(nevents,1); % for plot legend, etc
        frameseries = zeros(nevents,nframes);
        
        % extract timing data for events -- note TR conversion to frames
        
        for eindex = 1:numel(event_list)
            
            this_event = event_list(eindex);
            
            event_names{eindex} = this_event.name;           
            event_onsets = round(this_event.ons/TR);
            event_durations = round(this_event.dur/TR);      
           
            % fill frameseries by extending onsets by durations
            
            temp = zeros(1,nframes);
             
            for oindex = 1:numel(event_onsets)
                istart = event_onsets(oindex);
                if (istart==0);istart=1;end
                iend = istart + event_durations(oindex) - 1;
                if (iend<istart);iend=istart;end
                temp(istart:iend) = 1;
            end
            
            temp = temp(1:nframes); % in case last event ran over...
            
            frameseries(eindex,:) = temp;

        end
        
        % guess: we can return now -- there is only entry per subj/sess
            
        return;
 
    end
    
end

end
            
            
            