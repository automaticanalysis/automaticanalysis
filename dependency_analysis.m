clear; clc; % Clear the screen before proceeding

% We need to be inside the toolbox to work on it
cd(fileparts(which('aa_ver5')))

toolboxPath = pwd;

fldrDir = genpath(toolboxPath);
addpath(fldrDir); % To add the path to this toolbox!

ind = 0;

Dependencies = [];
Dependencies.name = [];
Dependencies.deps = [];
Dependencies.depsI = [];
Dependencies.depsT = [];
Dependencies.number = [];

% Then recurse inside each directory until you run out of paths
while ~isempty(strtok(fldrDir, ':'))
    % Get each of the directories made by gendir
    [fldrCurr fldrDir] = strtok(fldrDir, ':');
    
    % Get all .m files in this folder
    
    D = dir(fullfile(fldrCurr, '*.m'));
    for d = 1:length(D)
        if ~strcmp(D(d).name, [mfilename '.m'])
            ind = ind + 1;
            
            Dependencies(ind).name = fullfile(fldrCurr, D(d).name);
        end
    end
end

%% Find out MATLAB Toolbox dependecies
h = waitbar(0,'Cheked for MATLAB Toolbox dependecies');
for ind = 1:numel(Dependencies)    
    clear deps
    try
        deps = dependencies.toolboxDependencyAnalysis({Dependencies(ind).name});
    catch
        deps = 'failed';
    end
    
    deps(strcmp(deps,'MATLAB')) = [];
    deps(strcmp(deps,'Statistical Parametric Mapping')) = [];
    deps(strcmp(deps,'TFCE Toolbox')) = [];    
    
    Dependencies(ind).depsT = deps;
    waitbar(ind/numel(Dependencies),h);
end
close(h);

%% Find out dependency!
h = waitbar(0,'Cheked for file dependecies');
for ind = 1:length(Dependencies)    
    fprintf('\nWorking %d/%d: %s', ind, length(Dependencies), Dependencies(ind).name)
    
    clear deps
    deps = depfun(Dependencies(ind).name, '-quiet');
    pause(0.0001);
    
    for e = length(deps):-1:1
        % If toolbox path is not in the dependency, ignore it!
        if isempty(strfind(deps{e}, toolboxPath))
            deps(e) = [];
        elseif strcmp(deps{e}, Dependencies(ind).name)
            deps(e) = [];
        end
    end
    Dependencies(ind).deps = deps;
    waitbar(ind/numel(Dependencies),h);
end
close(h);

%% Find out inverse dependency!
h = waitbar(0,'Cheked for inverse dependecies');
for ind = 1:length(Dependencies)
    Dependencies(ind).number = 0;
    Dependencies(ind).depsI = {};
    
    for d = 1:length(Dependencies)
        if d ~= ind
            for e = 1:length(Dependencies(d).deps)
                if strcmp(Dependencies(ind).name, Dependencies(d).deps{e})
                    Dependencies(ind).depsI(end+1) = {Dependencies(d).name};
                    Dependencies(ind).number = Dependencies(ind).number + 1;
                end
            end
        end
    end
    waitbar(ind/numel(Dependencies),h);
end
close(h);

unused_deps = {};
fprintf('\n')

%% Find out any functions that are not used by anything!
h = waitbar(0,'Cheked for orphan functions');
for ind = 1:length(Dependencies)
    [path name] = fileparts(Dependencies(ind).name);
    if Dependencies(ind).number == 0 && ... % The number of uses of this function must be 0 in total...
            isempty(strfind(name, 'aamod')) && ... % ...and it must not be a module...
        isempty(strfind(name, 'aa_user')) % ... nor a user script
        unused_deps{end+1} = Dependencies(ind).name;
        
        fprintf([Dependencies(ind).name '\n']);
    end         
    waitbar(ind/numel(Dependencies),h);
end
close(h);
%% Now find how many dependencies each file has...
save('Dependencies.mat', 'Dependencies', 'unused_deps');