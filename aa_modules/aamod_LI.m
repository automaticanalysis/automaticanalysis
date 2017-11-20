function [aap,resp] = aamod_LI(aap, task, subj, sess)
%
% aamod_LI
%
% compute laterality index using Markos Wilke's LI toolbox
%
% See: 'A new toolbox to assess lateralization in functional MR-data', J Neurosci Meth, 2007, 163:128-136
%
% To use, download the toolbox from (as of this writing):
%
%	https://www.medizin.uni-tuebingen.de/kinder/en/research/neuroimaging/software/?download=li-toolbox
%
% This will unzip as a folder named LI_toolbox. Move this folder to /path/to/your/spm/toolbox and
% restart aa (the LI toolbox will also appear from the Toolbox pull-down menu in the SPM GUI, if
% you'd like to use it that way instead).
%
% See LI_manual.pdf included in the download for command options. This pdf also includes preprints of
% articles describing the algorithm. Note Wilke requests developers do not package his toolkit in other
% software, which is why it's not included in the aa distribution.
%

% init

resp='';

switch task
	
    case 'domain'
		
        resp='session'; 
		
    case 'report'
		
		% add txt results; add graphics results if jpeg versions exist in the directory
		
		temp_string = 'Lateralization Statistics';
		
		savedir = pwd;
		session_path = aas_getsesspath(aap, subj, sess);
		cd(session_path);
		
		% the default filename 'li.txt' is hardcoded in LI toolbox
				
		lstats_fname = fullfile(session_path, 'li.txt');
		
		fid = fopen(lstats_fname,'r');
		
		if (fid < 0)
			aap = aas_report_add(aap, '<h4>Cannot find lateralization analysis summary.</h4>');
			return;
		end
		
		while (ischar(temp_string))
			aap = aas_report_add(aap, subj, sprintf('<h3>%s</h3>', temp_string));
			temp_string = fgetl(fid); % will return -1 at EOF
		end
		
		fclose(fid);
				
		% LI toolbox creates postscript graphics by default. These can't be embedded in the report because HTML, which 
		% leaves two options: 1) convert the ps to jpeg here which requires, e.g., GhostScript because Matlab doesn't
		% doesn't provide a utility for it (thanks Mathworks), or 2) edit your installation of LI_toolbox to create
		% jpeg instead of postscript. The latter is easily done by changing the switch in the print commands in LI.m, 
		% LI_iter.m, and LI_boot.m, from -dpsc2 to -djpeg. But if we don't find any jpegs, and gs conversion fails,
		% print a message to help the user to get their ducks aligned.
		
		jpg_files = dir(fullfile(session_path,'*.jpg'));
		ps_files = dir(fullfile(session_path,'*.ps'));
					
		if (numel(jpg_files) == 0 && numel(ps_files) > 0)
			
			% assume the go-to converter is ghostscript
			
			[ ierr, gs ] = system('which gs');
			
			% sometimes unix $PATH doesn't get loaded into matlab correctly
			% -- if we didn't find gs, try looking for it again where it probably is
			
			if (ierr)
				[ ierr, gs ] = system('which /usr/local/bin/gs');
			end
			
			if (ierr)
				
				% oh well, we tried...
			
				aas_log(aap, false, '\nCannot add LI graphics to the report because they are in postscript.');
				aas_log(aap, false, 'If you want to have LI_toolbox graphics included to the report, change the print');
				aas_log(aap, false, 'statements in LI.m, LI_iter.m, and LI_boot.m to use -djpeg instead of -dpsc2');
				aas_log(aap, false, 'or install a postscript-to-jpeg converter such as Ghostscript.\n');
				return;
				
			else
				
				aas_log(aap, false, '\nConverting postscript output to jpeg to include in report...\n');
				
				for index = 1:numel(ps_files)
					
					[ p,name,e ] = fileparts(ps_files(index).name);

					% gs conversion gives better results if you tell it the dimensions of the image
					
					gs = deblank(gs);
					
					command = sprintf('%s -sDEVICE=bbox -dNOPAUSE -dBATCH -q %s.ps', gs, name);
					[ ierr, result ] = system(command);
					
					if (~ierr)
						key = 'BoundingBox';
						index = min(strfind(result,key));
						temp = sscanf(result(index+length(key)+1:end),'%g',4);
						w = round(1.2 * (double(temp(1)) + double(temp(3))));
						h = round(1.2 * (double(temp(2)) + double(temp(4))));
						if (isnan(w) || isnan(h)) % sanity check
							command = sprintf('%s -sDEVICE=jpeg -dJPEGQ=100 -dNOPAUSE -dBATCH -dSAFER -q -sOutputFile=%s.jpg %s.ps', gs, name, name);
						else
							command = sprintf('%s -sDEVICE=jpeg -dJPEGQ=100 -dNOPAUSE -dBATCH -g%dx%d -dSAFER -q -sOutputFile=%s.jpg %s.ps', gs, w, h, name, name);
						end
						system(command);
					end
					
					
				end
				
			end

		end
		
		% verify jpeg list (may have been updated by ps conversion)
		
		jpg_files = dir(fullfile(session_path,'*.jpg'));
			
		for index = 1:numel(jpg_files)
			aap = aas_report_add(aap, subj, '<table><tr><td>');
			aap = aas_report_addimage(aap, subj, fullfile(session_path, jpg_files(index).name));
			aap = aas_report_add(aap, subj,'</td></tr></table>');
		end			

		% pop the savedir
		
		cd(savedir);
		
    case 'doit'
	
		% the LI toolbox creates a bunch of files in the working dir
		% we don't have any control over their code so cd to the sess dir
		
		savedir = pwd;
		session_path = aas_getsesspath(aap, subj, sess);
		cd(session_path);
		
		% extract the options into an LI struct
		
		LI_command_struct.B1 = aap.tasklist.currenttask.settings.LI_inclusive_mask;
		LI_command_struct.C1 = aap.tasklist.currenttask.settings.LI_exclusive_mask;
		LI_command_struct.thr1 = aap.tasklist.currenttask.settings.LI_thresholding;
	
		% the only thing missing from the command string is the file(s) to process
		
		input_streamname = aas_getstreams(aap,'input');
		fname_list = aas_getfiles_bystream(aap, subj, sess, input_streamname{1});
	
		for index = 1:size(fname_list,1)
			LI_command_struct.A = fname_list(index,:);
			LI(LI_command_struct);
		end
		
		% the LI toobox leaves the SPM graphics window up. Take it down.
	
		spm_figure('Close', 'Graphics'); 
		
        %  describe outputs
				
		% LI creates (at least) three output files (more, depending on options):
		%
		%	LI_boot.ps
		%	LI_masking.ps
		%	li.txt
		%
		% we'll save the stats as an output stream. The others
		% we just save to the report (if they've been jpeg converted)					
		
 		LI_stats_fname = fullfile(session_path, 'li.txt');
		aap = aas_desc_outputs(aap, subj, sess, 'LI_stats', LI_stats_fname);
				
		% pop the sess dir
		
		cd(savedir);
		
		
    case 'checkrequirements'
					
        
    otherwise
        aas_log(aap, 1, sprintf('%s: Unknown task %s', mfilename, task));
		
		
end

