function search_code(snippet, editFile)
if nargin < 2
    editFile = 0;
end

% We need to be inside the toolbox to work on it
cd(fileparts(mfilename('fullpath')))

toolboxPath = pwd;

fldrDir = genpath(toolboxPath);
addpath(fldrDir); % To add the path to this toolbox!

ind = 0;


fprintf('\nSnippet %s is found within...\n', snippet)
% Then recurse inside each directory until you run out of paths
while ~isempty(strtok(fldrDir, ':'))
    % Get each of the directories made by gendir
    [fldrCurr fldrDir] = strtok(fldrDir, ':');
    
    % Get all .m files in this folder
    D = dir(fullfile(fldrCurr, '*.m'));
    for d = 1:length(D)
        T = textread(fullfile(fldrCurr, D(d).name), '%s', 'whitespace', '', 'bufsize', 1024^2);
        
        if ~isempty(strfind(T{1}, snippet))
           fprintf('%s\n', fullfile(fldrCurr, D(d).name)) 
           if editFile
               edit(fullfile(fldrCurr, D(d).name))
           end
        end
    end
    % And all .xml files
    D = dir(fullfile(fldrCurr, '*.xml'));
    for d = 1:length(D)
        T = textread(fullfile(fldrCurr, D(d).name), '%s', 'whitespace', '', 'bufsize', 1024^2);
        
        if ~isempty(strfind(T{1}, snippet))
           fprintf('%s\n', fullfile(fldrCurr, D(d).name))
           if editFile
               edit(fullfile(fldrCurr, D(d).name))
           end
        end
    end
end
