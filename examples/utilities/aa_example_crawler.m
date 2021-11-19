function aa_example_crawler(keyword)

% this function will crawl $AAHOME/example and echo the header 
% block of all mfiles to the command window. Optionally, only the
% mfiles containing the search string "keyword" will be displayed

if ~exist('keyword','var') || isempty(keyword)
    keyword = '';
end

% get a list of aa example scripts

[ aaInstallDir,~,~ ] = fileparts(which('aa_ver5'));

exampledir = fullfile(aaInstallDir,'examples');
command = sprintf('find %s -name \\*.m', exampledir);

[ status,script_list ] = system(command);

if (status || ~numel(script_list))
	error('script list generation failed');
end

script_list = split(deblank(script_list));

moreState = get(0,'More');
more on

for index = 1:numel(script_list)
    
    if ~isempty(keyword)    
        command = sprintf('grep %s %s', keyword, script_list{index});        
        [ status,~ ] = system(command);
        if (status > 0)
            % grep didn't find search string in the file -- don't display:
            continue;
        end
    end
    
    [~,scriptname,~] = fileparts(script_list{index});
    fprintf('\n <strong> Script: %s </strong> \n\n', scriptname);
    help(script_list{index})
end

% don't turn more off is user had it on
if strcmp(moreState,'off');more off;end


