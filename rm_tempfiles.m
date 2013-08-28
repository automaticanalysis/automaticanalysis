function rm_tempfiles
% Again we need to be inside the toolbox to work on it
cd(fileparts(mfilename('fullpath')))

try
    % First get the path inside the folders here
    fldrDir = genpath(pwd);
    % We don't want stuff ending in...
    DONTwant = {'.m~', '.xml~', 'txt~'};
    % Then recurse inside each directory until you run out of paths
    while ~isempty(strtok(fldrDir, ':'))
        % Get each of the directories made by gendir
        [fldrCurr fldrDir] = strtok(fldrDir, ':');
        
        % Remove any temporary files we might have created
        if ~isempty(dir([fldrCurr '/*.m~']))
            unix(['rm -rf ' fldrCurr '/*.m~']);
        end
        if ~isempty(dir([fldrCurr '/*.xml~']))
            unix(['rm -rf ' fldrCurr '/*.xml~']);
        end
        if ~isempty(dir([fldrCurr '/*.txt~']))
            unix(['rm -rf ' fldrCurr '/*.txt~']);
        end
    end
end

end
