function aa_jpeg_crawler(aap, path_regex)
%
% crawl an aa results directory for jpegs and assemble them into
% a scrollable html doc
%
%	usage: aa_jpeg_crawler(aap, [path_regex])
%
% INPUT
%
%	aap - aap struct of analysis to crawl
%
% OPTIONAL
%
%   path_regex - only include files in the report having this regex
%   in their pathname (e.g. "firstlevel_model" will only include files
%   *firstlevel_model*.jpg in the doc)
%
% Usage
%
% aa_jpeg_crawler crawls an analysis directory and assembles all jpgs
% into a single scrollable html document for easy review. A clickable
% table-of-contents is generated in the left column that will jump to
% results from the selected module. A preamble section includes the
% provenance map (aap_prov.dot). Not dot must be installed somewhere in
% the shell path to convert the provenance map to jpeg for display,
% otherwise it will be skipped.
%
% Clicking on any image in the report will open the it in a new tab
% window, which can be zoomed for inspection. The image path is displayed
% as the browswer address bar, as well appearing as a tooltip when the
% mouse is hovered over the image.
%
% CHANGE HISTORY
%
% 12/2019 - change exclude list to include (path_regex);
%           add "no images found" to empty section headers
% 11/2019 - include both .jpg and .jpeg files
% 03/2019 - add exclusion list
% 07/2018 - make sections collapsible
% 05/2018 - new
%

savedir = pwd;

if (isempty(aap))
	disp('Usage: aa_jpeg_crawler(aap,[path_regex])');
	return;
end

if (nargin < 2)
    path_regex = '';
end

cd([aap.acq_details.root '/' aap.directory_conventions.analysisid]);

% first thing, see if provenance exists and needs converted
% UPDATE: best to delete the old jpg so we force it to be current

system('rm -f aap_prov.jpg');

if (~exist('aap_prov.jpg','file'))
% 	system('dot -Tjpg aap_prov.dot > aap_prov.jpg 2>/dev/null');
	system('dot aap_prov.dot -Tjpg -o aap_prov.jpg 2>/dev/null');
end

% saving throw: if aap_prov.jpg still doesn't exist, maybe Matlab
% can't find dot. Try looking one other place for it...

if (~exist('aap_prov.jpg','file'))
	system('/usr/local/bin/dot aap_prov.dot -Tjpg -o aap_prov.jpg 2>/dev/null');
end

% make a list of all the jpg/jpegs under the directory
%
% sorting on field #9 works well, i.e.
%
% status = system('find -s `pwd` -name \*.jpg | sort -t/ -k9 > montage_temp.txt');
%
% (this is actually root depth + 4)

root_depth = numel(split(aap.acq_details.root,'/'));
    
command = sprintf('find -s -E `pwd` -type f  -path ''*%s*''  -iregex ''.*[Jj][Pp][Ee]?[Gg]'' | sort -t/ -k%d > montage_temp.txt', path_regex, root_depth+4);

status = system(command);

if (status); cleanup_and_exit(1); end

fid = fopen('aa_jpeg_crawler.htm','w');
if (fid < 0); cleanup_and_exit(2); end

html_make_head(fid);
html_add_table_of_contents(fid, aap);
html_open_div(fid);
html_add_jpegs(fid, aap);
html_close_div(fid);
html_close(fid);

cleanup_and_exit(0);
cd(savedir);

%-----------------------------------------------------------------------------------------------------------------------------------
% cleanup_and_exit - THIS MUST BE NESTED FOR VARIABLE ACCESS
%-----------------------------------------------------------------------------------------------------------------------------------

function cleanup_and_exit(ierr)

	system('rm -f montage_temp01.txt');
	system('rm -f montage_temp.txt');
	if exist('fid');fclose(fid);end
	cd(savedir);
	if (ierr)
		if exist('aap','var')
			aas_log(aap, true, sprintf('\n%s: Montage generation failed (ierr = %d.\n', mfilename, ierr));
		else
			error('Montage generation failed (ierr = %d).\n', ierr);
		end
	end
	
end


end % NESTED


%-----------------------------------------------------------------------------------------------------------------------------------
% html_make_head
%-----------------------------------------------------------------------------------------------------------------------------------

function html_make_head(fid)

	% create the html head, which now includes css for index
	
	% <meta name="viewport" content="width=device-width, initial-scale=1">

	fprintf(fid,'%s\n','<!DOCTYPE html>');
	fprintf(fid,'%s\n','<html>');
	fprintf(fid,'%s\n','<head>');
	fprintf(fid,'%s\n','<style>');
	fprintf(fid,'%s\n','body {');
	fprintf(fid,'%s\n','    margin: 0;');
	fprintf(fid,'%s\n','}');
	fprintf(fid,'%s\n','');
	fprintf(fid,'%s\n','ul {');
	fprintf(fid,'%s\n','    list-style-type: none;');
	fprintf(fid,'%s\n','    margin: 0;');
	fprintf(fid,'%s\n','    padding: 0;');
	fprintf(fid,'%s\n','    width: 12%;');
	fprintf(fid,'%s\n','    background-color: #f1f1f1;');
	fprintf(fid,'%s\n','    position: fixed;');
	fprintf(fid,'%s\n','    height: 100%;');
	fprintf(fid,'%s\n','    overflow: auto;');
	fprintf(fid,'%s\n','}');
	fprintf(fid,'%s\n','');
	fprintf(fid,'%s\n','li a {');
	fprintf(fid,'%s\n','    display: block;');
	fprintf(fid,'%s\n','    color: #000;');
	fprintf(fid,'%s\n','    padding: 8px 16px;');
	fprintf(fid,'%s\n','    text-decoration: none;');
	fprintf(fid,'%s\n','    font-size:12px;');
	fprintf(fid,'%s\n','    font-weight:bold;');
	fprintf(fid,'%s\n','}');
	fprintf(fid,'%s\n','');
	fprintf(fid,'%s\n','li a.active {');
	fprintf(fid,'%s\n','    background-color: #4CAF50;');
	fprintf(fid,'%s\n','    color: white;');
	fprintf(fid,'%s\n','}');
	fprintf(fid,'%s\n','');
	fprintf(fid,'%s\n','li a:hover:not(.active) {');
	fprintf(fid,'%s\n','    background-color: #555;');
	fprintf(fid,'%s\n','    color: white;');
	fprintf(fid,'%s\n','}');
	
    fprintf(fid,'%s\n','.collapsible {');
    fprintf(fid,'%s\n','    background-color: #777;');
    fprintf(fid,'%s\n','    color: white;');
    fprintf(fid,'%s\n','    cursor: pointer;');
    fprintf(fid,'%s\n','    padding: 18px;');
    fprintf(fid,'%s\n','    width: 100%;');
    fprintf(fid,'%s\n','    border: none;');
    fprintf(fid,'%s\n','    text-align: left;');
    fprintf(fid,'%s\n','    outline: none;');
    fprintf(fid,'%s\n','    font-size: 15px;');
    fprintf(fid,'%s\n','}');
    fprintf(fid,'%s\n','');
    fprintf(fid,'%s\n','.active, .collapsible:hover {');
    fprintf(fid,'%s\n','    background-color: #555;');
    fprintf(fid,'%s\n','}');
    fprintf(fid,'%s\n','');
    fprintf(fid,'%s\n','.content {');
    fprintf(fid,'%s\n','    padding: 0 18px;');
    fprintf(fid,'%s\n','    display: block;');
    fprintf(fid,'%s\n','    overflow: hidden;');
    fprintf(fid,'%s\n','    background-color: #f1f1f1;');
    fprintf(fid,'%s\n','}');
	
	fprintf(fid,'%s\n','</style>');
	fprintf(fid,'%s\n','</head>');
	fprintf(fid,'%s\n','<body>');

	fprintf(fid,'%s\n','<title>aa Report</title>');

end

%-----------------------------------------------------------------------------------------------------------------------------------
% html_add_table_of_contents
%-----------------------------------------------------------------------------------------------------------------------------------

function html_add_table_of_contents(fid, aap)

	% add a clickable table of contents tab based on modules
	
	% the entries go like:
	%
	% 	<li><a href="#tag01">aamod_firstlevel_threshold_00001</a></li>
	% 	<li><a href="#tag02">aamod_norm_noss_00001</a></li>
	% 	<li><a href="#tag03">aamod_firstlevel_model_00001</a></li>
	% 	<li><a href="#tag03">aamod_firstlevel_model_00002</a></li>
	% 	<li><a href="#tag03">aamod_firstlevel_model_00003</a></li>

	% call this after make_html_head 

	% this is based on aap if module_list is a struct
	% or parse out the modules if module_list is a cell array

	fprintf(fid,'<ul>\n');
	fprintf(fid,'<br/> &nbsp; CONTENTS<br/><br/>\n');

	% first entry is provenance if a jpeg of the map exists
		
	if (exist('aap_prov.jpg','file'))
		fprintf(fid,'<li><a href="#provenance">provenance</a></li>\n');
	end
	
	% also include any diagnostic jpegs present at the root level
	
	root_diagnostic_files = dir('diagnostics_*.jpg');

	for index = 1:numel(root_diagnostic_files)
		fname = root_diagnostic_files(index).name;
		[~,label,~] = fileparts(fname);
		fprintf(fid,'<li><a href="#%s">%s</a></li>\n', label, label);
	end
	
	fprintf(fid,'<hr/>\n');
	 
	tasklist = aap.tasksettings;
	tasknames = fieldnames(tasklist);
	tagnum = 1;
	
	for index = 1:numel(tasknames)
		
		% handle possibility of multiple module appearance
		
		module_name = tasknames{index};
		module_count = numel(tasklist.(module_name));
		
		for count = 1:module_count
			fprintf(fid,'<li><a href="#tag%02d">%s_%05d</a></li>\n', tagnum, module_name, count);
			tagnum = tagnum + 1;
		end
		
	end
				
  	fprintf(fid,'</ul>\n');

end


%-----------------------------------------------------------------------------------------------------------------------------------
% html_open_div
%-----------------------------------------------------------------------------------------------------------------------------------

function html_open_div(fid)

	% need an extra div tag here to balance the first call to
	% collapsible button...
	
	fprintf(fid, '<div style="margin-left:12%%;padding:1px 16px;height:1000px;"><div>\n');

end


%-----------------------------------------------------------------------------------------------------------------------------------
% html_add_jpegs
%-----------------------------------------------------------------------------------------------------------------------------------

function html_add_jpegs(fid, aap)

	if (exist('aap_prov.jpg','file'))
		fprintf(fid,'<h2 id="provenance">provenance</h2>\n'); 
		html_add_image(fid, fullfile(pwd,'aap_prov.jpg'));
	end
	
	% also include any diagnostic jpegs present at the root level
	
	root_diagnostic_files = dir('diagnostics_*.jpg');

	for index = 1:numel(root_diagnostic_files)
		fname = root_diagnostic_files(index).name;
		[~,label,~] = fileparts(fname);
		fprintf(fid,'<h2 id="%s">%s</h2>\n',label,label); 
		html_add_image(fid, fullfile(pwd,fname));
	end

	% crawl tasklist, add jpegs

	tasklist = aap.tasksettings;
	tasknames = fieldnames(tasklist);
	tagnum = 1;

	for index = 1:numel(tasknames)
		
		current_image_type = [];

		% handle multiple module appearance
		
		module_name = tasknames{index};
		module_count = numel(tasklist.(module_name));
		
		for count = 1:module_count
			
			stagename = sprintf('%s_%05d', module_name, count);
			
			[ status, jpeg_list ] = system(sprintf('grep %s montage_temp.txt', stagename));

			% grep returns status 1 if nothing found (but its not an error)

			if (status>1); cleanup_and_exit(3); end
            
            if (isempty(jpeg_list))
                fprintf(fid,'<h2 id="tag%02d">%s (no images found)</h2>\n', tagnum, stagename);
            else
                fprintf(fid,'<h2 id="tag%02d">%s</h2>\n', tagnum, stagename);
            end
            
			tagnum = tagnum + 1;

			%  split adds a (wrong) entry if jpeg_list is empty...

			if (~isempty(jpeg_list));jpeg_list = split(deblank(jpeg_list));end

			% jpeg_list should be sorted by subject. We put all images of a type in a
			% row and insert a linebreak each time we encounter a new type
			% TEST ME -- this may not be needed with tasklist-based sorting
			
			for jindex = 1:numel(jpeg_list)

				jpeg_name = jpeg_list{jindex};

				[ ~,image_type,~ ] = fileparts(jpeg_name);

				if (isempty(current_image_type))
					current_filetype = image_type;
				end
				
				if ~strcmp(image_type, current_image_type)
					html_add_linebreak(fid,image_type);
					current_image_type = image_type;
				end

				html_add_image(fid, jpeg_name);

			end
			
			fprintf('.');

		end

	end
	
	fprintf('\n');
	
end


%-----------------------------------------------------------------------------------------------------------------------------------
% html_add_image
%-----------------------------------------------------------------------------------------------------------------------------------

function html_add_image(fid, image_path)

	% linux is doing this weird thing where it creates corrupt copies
	% of jpegs we create ; the names begin with "._" - for now,
	% just bounce these.
	
	if strfind(image_path,'/.')
		return;
	end

	h = imfinfo(image_path);

	% REVISIT: scale factor by numel(jpeg_list)?

	scaling_factor = 500;

	if (h.Width > h.Height)
		width = scaling_factor;
		height = round(h.Height * scaling_factor / h.Width);
	else
		height = scaling_factor;
		width = round(h.Width * scaling_factor / h.Height);
	end

	fprintf(fid, sprintf('<a href="file://%s" target="_blank">\n', image_path));
	fprintf(fid, sprintf('<img src="%s" title="%s" height="%d" width="%d"></a>\n', image_path, image_path, height, width));

end
			

%-----------------------------------------------------------------------------------------------------------------------------------
% html_add_linebreak
%-----------------------------------------------------------------------------------------------------------------------------------

function html_add_linebreak(fid,s)

	if (isempty(s))
		fprintf(fid, '<br/>\n');
	else
		fprintf(fid, '</div><br/><br/><button class="collapsible">%s</button>\n', s);
		fprintf(fid, '<div class="content"><br/>\n');
    end
	
end

%-----------------------------------------------------------------------------------------------------------------------------------
% html_close_div
%-----------------------------------------------------------------------------------------------------------------------------------

function html_close_div(fid)
	fprintf(fid, '</div><br/><br/>');
end

%-----------------------------------------------------------------------------------------------------------------------------------
% html_close
%-----------------------------------------------------------------------------------------------------------------------------------

function html_close(fid)


fprintf(fid,'%s\n','<script>');
fprintf(fid,'%s\n','var coll = document.getElementsByClassName("collapsible");');
fprintf(fid,'%s\n','var i;');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','for (i = 0; i < coll.length; i++) {');
fprintf(fid,'%s\n','    coll[i].addEventListener("click", function() {');
fprintf(fid,'%s\n','        this.classList.toggle("active");');
fprintf(fid,'%s\n','        var content = this.nextElementSibling;');
fprintf(fid,'%s\n','        if (content.style.display === "none") {');
fprintf(fid,'%s\n','            content.style.display = "block";');
fprintf(fid,'%s\n','        } else {');
fprintf(fid,'%s\n','            content.style.display = "none";');
fprintf(fid,'%s\n','        }');
fprintf(fid,'%s\n','    });');
fprintf(fid,'%s\n','}');
fprintf(fid,'%s\n','</script>');

fprintf(fid, '</body></html>');

end




