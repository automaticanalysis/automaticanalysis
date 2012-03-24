function commitSVN(SVNmssg)
% Again we need to be inside the toolbox to work on it
cd(fileparts(mfilename('fullpath')))

try
    % Let's add all the files present in the toolbox
    % First get the path inside the folders here
    fldrDir = genpath(pwd);
    % We don't want stuff ending in...
    DONTwant = {'.m~', '.xml~', 'txt~'};
    % Then recurse inside each directory until you run out of paths
    while ~isempty(strtok(fldrDir, ':'))
        % Get each of the directories made by gendir
        [fldrCurr fldrDir] = strtok(fldrDir, ':');
        % Check it's not a .svn folder
        if isempty(strfind(fldrCurr, '.svn'))
            D = dir(fldrCurr);
            for d = 1:length(D)
                if ~D(d).isdir || ~any(strfind(D(d).name, '.'))
                    [~, statusS] = unix(['svn status ' fldrCurr '/' D(d).name]);
                    if ~isempty(statusS)
                        fprintf('\n%s', statusS);
                        if strcmp(statusS(1), 'M')
                            % Modified file, don't do anything
                        elseif strcmp(statusS(1), 'D')
                            % Deleted file, don't do anything
                        elseif strcmp(statusS(1), 'A')
                            % Already added
                        elseif strcmp(statusS(1), '?')
                            % Don't know what this is, let's try to add it
                            % Add files files present here
                            addFile(fldrCurr, D(d).name, DONTwant)
                        end
                    else
                        % Unmodified file, don't do anything
                    end
                else
                    % It is one of the . or .. folders
                end
            end
            
            % Remove any temporary files we might have created
            if ~isempty(dir([fldrCurr '/*.m~']))
                unix(['svn rm --force ' fldrCurr '/*.m~']);
            end
            if ~isempty(dir([fldrCurr '/*.xml~']))
                unix(['svn rm --force ' fldrCurr '/*.xml~']);
            end
            if ~isempty(dir([fldrCurr '/*.txt~']))
                unix(['svn rm --force ' fldrCurr '/*.txt~']);
            end
        end
    end
    
    % Finally, all we need to do, is to commit the changes
    unix(['svn commit -m ''' SVNmssg '''']);
catch
    error('If this does not work you either...\n1) are in the wrong directory\n2) have not added an update message')
end
end

function addFile(fldrCurr, Fname, DONTwant)
want = 1;
for w = 1:length(DONTwant)
    if length(Fname)>=length(DONTwant{w})
        if strcmp(DONTwant{w}, ...
                Fname((end+1-length(DONTwant{w})):end))
            want = 0;
        end
    end
end
if want
    unix(['svn add ' fldrCurr '/' Fname]);
end
end