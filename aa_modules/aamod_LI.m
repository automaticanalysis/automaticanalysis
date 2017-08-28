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
		
		session_path = aas_getsesspath(aap, subj, sess);
		
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
		
		% LI toolbox creates postscript graphics. These can't be embedded in the report because HTML, which leaves 
		% us two options: 1) convert the ps to jpeg here which requires, e.g., GhostScript because Matlab doesn't
		% doesn't provide a utility for it (thanks Mathworks), or 2) edit your installation of LI_toolbox to create
		% jpeg instead of postscript. The latter is easily done by changing the switch in the print commands in LI.m, 
		% LI_iter.m, and LI_boot.m, from -dpsc2 to -djpeg. Let's assume if the user really wants LI output in the
		% report, they're savvy enough to do a simple code edit.
	
		% so here we just add every jpg in the session folder to the report (these can vary depending on which options
		% were used). If there are no jpeg, we add a message to help the user to get their ducks lined up	
		
		jpeg_files_present = dir(fullfile(session_path,'*.jpg'));
		
        if (numel(jpeg_files_present) == 0)
			aap = aas_report_add(aap, subj, '<h3>Cannot add LI graphics to the report because they are in postscript.</h3>');
			aap = aas_report_add(aap, subj, '<h3>If you want to have LI_toolbox results included to the report, change the print</h3>');
			aap = aas_report_add(aap, subj, '<h3>statements in LI.m, LI_iter.m, and LI_boot.m to use -djpeg instead of -dpsc2</h3>');
			return;
		end
	
        for index = 1:numel(jpeg_files_present)
            aap = aas_report_add(aap, subj, '<table><tr><td>');
            aap = aas_report_addimage(aap, subj, fullfile(session_path, jpeg_files_present(index).name));
            aap = aas_report_add(aap, subj,'</td></tr></table>');
		end
		
         
    case 'doit'
	
		% the LI toolbox creates a bunch of files in the working dir
		% we don't have any control over their code so cd to the sess dir
		
		savedir = pwd;
		session_path = aas_getsesspath(aap, subj, sess);
		cd(session_path);
		
		% LI has an odd API consisting of a struct built from lots of
		% required and optional command options. We just let the user
		% pass in a command string rather than individual params and
		% parse the string into a struct
		
		LI_command_string = aap.tasklist.currenttask.settings.LI_command_string;		
		temp = strsplit(LI_command_string,',');
		fieldnames = { temp{1:2:numel(temp)} };
		fieldvals =  { temp{2:2:numel(temp)} };
		fieldvals =  num2cell(cellfun(@str2num,fieldvals));
		LI_command_struct = cell2struct(fieldvals,fieldnames,2);

		% the only thing missing from the command string is the file to process
		
		input_streamname = aas_getstreams(aap,'input');
		fname_list = aas_getfiles_bystream(aap, subj, input_streamname{1});
	
		for index = 1:size(fname_list,1)
			LI_command_struct.A = fname_list(index,:);
			LI(LI_command_struct);
		end
		
        %  describe outputs
				
		% LI creates (at least) three output files (more, depending on options):
		%
		%	LI_boot.ps
		%	LI_masking.ps
		%	li.txt
		%
		% we'll save the stats (li.txt) as an output stream. The others
		% we just save to the report (if they've been jpeg converted)					
		
 		LI_stats_fname = fullfile(session_path, 'li.txt');
		aap = aas_desc_outputs(aap, subj, sess, 'LI_stats', LI_stats_fname);
				
		% pop the sess dir
		
		cd(savedir);
		
		
    case 'checkrequirements'
			
		if isempty(which('LI'))
			aas_log(aap, true, sprintf('\n%s: LI not found. Make sure it exists in /path/to/spm/toolbox.\n', mfilename));
		end
		
        
    otherwise
        aas_log(aap, 1, sprintf('%s: Unknown task %s', mfilename, task));
		
		
end