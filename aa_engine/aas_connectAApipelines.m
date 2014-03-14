function aap = aas_connectAApipelines(aap, remoteAAlocations)
% aas_connect_remoteAA(aap, remoteAAlocations)
% Connects streams from other (remote) AA pipelines to the local one.
%
% The "local" AA analysis is the one that is currently running (aap input)
%
% The "remote" AA analysis lives somewhere else (the remoteAAPfn input).
% If this is a cell array of diretories (i.e., multiple AA locations), then
% priority is given to AA locations earlier in the list.  If a stream isn't
% found in an earlier location, then it drops down to the next location.
%
% Any modules of the local AA that have input streams that aren't outputs
% of any previous modules in the local AA wil be connected up to the remote
% AA.

% Error checking:
if isempty(aap.acq_details.subjects)
    aas_log(aap, 1, 'aas_connectAApipelines() should be used after you have added subjects in your user script.');
end

if isempty(aap.acq_details.sessions)
    aas_log(aap, 1, 'aas_connectAApipelines() should be used after you have added sessions in your user script.');
end

if any(~isfield(remoteAAlocations, {'host', 'directory', 'allowcache', 'maxstagetag'}))
    aas_log(aap, 1, 'remoteAAlocations (input to aas_connectAApipelines) should be a struct array with the following fields: ''host'', ''directory'', ''allowcache'', ''maxstagetag''');
end

global aaworker;
try
    aaworker.parmpath;
catch
    [pth nme ext] = fileparts(tempname);
    aaworker.parmpath = aaworker_getparmpath(aap,[filesep sprintf('%s_%s',datestr(now,30),nme)]);
    aas_makedir(aap, aaworker.parmpath);
end

% We need to transfer over the remote AAP files, will put them here
studyPath = aas_getstudypath(aap);
if ~exist(studyPath), aas_makedir(aap, studyPath); end

% Collect remote AA structures here
remoteAA = {};

% Names, modules, and location indices of remote outputstreams
remoteOutputs = struct('name', {{}}, 'locI', {}, 'modI', {});

% Collect output streams from the remote AA locations
% Reverse order, because priority is given to earlier array elements in remoteAAlocations.
for locI = length(remoteAAlocations) : -1 : 1
    
    % There should be this file at the remote AA location
    remoteAAPfn = fullfile(remoteAAlocations(locI).directory, 'aap_parameters.mat');
    
    % This is where we will put it
    destAAPfn = fullfile(studyPath, sprintf('aap_parameters_%05d.mat', locI));
    
    % Copy the remote AAP file, could error here if the file can''t be found remotely
    aap = aas_copyfromremote(aap, remoteAAlocations(locI).host, remoteAAPfn, destAAPfn, 'allowcache', remoteAAlocations(locI).allowcache);
    
    % Load it!
    if ~exist(destAAPfn, 'file'), aas_log(aap, 1, sprintf('Cannot find %s', remoteAAPfn)); end
    remoteAA{locI} = load(destAAPfn); remoteAA{locI} = remoteAA{locI}.aap;
    
    % Check that all subjects in the local AA are present in the remote AA
    localSub = {aap.acq_details.subjects.mriname};
    remoteSub = {remoteAA{locI}.acq_details.subjects.mriname};
    subMatch = ismember(localSub, remoteSub);
    if any(~subMatch)
        aas_log(aap, 1, sprintf('Remote AAP (%s:%s) doesn''t have subjects: %s', remoteAAlocations(locI).host, remoteAAPfn, strjoin(localSub(~subMatch))));
    end
    
    % Check that all sessions in the local AA are present in the remote AA
    if ~(length(aap.acq_details.sessions)==1 && isempty(aap.acq_details.sessions(1).name))
    localSess = {aap.acq_details.sessions.name};
    remoteSess = {remoteAA{locI}.acq_details.sessions.name};
    sessMatch = ismember(localSess, remoteSess);
    if any(~sessMatch)
        aas_log(aap, 1, sprintf('Remote AAP (%s:%s) doesn''t have sessions: %s', remoteAAlocations(locI).host, remoteAAPfn, strjoin(localSess(~sessMatch))));
    end
    end
    
    maxModI = [];
    
    % Check for empty 'maxstagetag', otherwise get the module index from
    % the remoteAAP
    if isempty(remoteAAlocations(locI).maxstagetag)
        maxModI = length(remoteAA{locI}.tasklist.main.module);
    else
        [maxModI name id] = aas_getmoduleindexfromtag(remoteAA{locI}, remoteAAlocations(locI).maxstagetag);
    end
    
    % Find the last instance of all remote output streams
    for modI = 1 : maxModI
        mod = remoteAA{locI}.tasklist.main.module(modI);
        
        if isfield(remoteAA{locI}.tasksettings.(mod.name)(mod.index), 'outputstreams')
            outputStreams = remoteAA{locI}.tasksettings.(mod.name)(mod.index).outputstreams.stream;
            if ~iscell(outputStreams), outputStreams = {outputStreams}; end
            
            for oI = 1 : length(outputStreams)
                streamI = strcmp(outputStreams{oI}, {remoteOutputs.name});
                
                % If this output hasn't appeared yet, add it to the list
                if ~any(streamI)
                    remoteOutputs(end+1) = struct('name', outputStreams{oI}, 'locI', locI, 'modI', modI);
                    
                    % Otherwise, update the module and location indices
                else
                    remoteOutputs(streamI).modI = modI;
                    remoteOutputs(streamI).locI = locI;
                end
                
            end
        end
    end
end

% Track the names of output streams from modules in the local analysis.  If
% a stream is an output from a previous stage in the local AA, then we
% don't bother trying to connect it from the remote AA.
prevOutputs = {};

for modI = 1 : length(aap.tasklist.main.module)
    mod = aap.tasklist.main.module(modI);
    
    if isfield(aap.tasksettings.(mod.name)(mod.index).inputstreams, 'stream')
        
        % Names of input and output streams for this module
        inputStreams = aap.tasksettings.(mod.name)(mod.index).inputstreams.stream;
        if ~iscell(inputStreams), inputStreams = {inputStreams}; end
        
        remoteStreams = struct('stream', {}, 'stagetag', {}, 'sourcedomain', {}, 'host', {}, 'aapfilename', {}, 'allowcache', {});
        
        for iI = 1 : length(inputStreams)
            
            % If the input doesn't come from a previous module, let's add it to
            % the list of remote streams for this module.
            if ~ismember(inputStreams{iI}, prevOutputs)
                
                % Is the input stream present in the list of remote streams?
                rI = strcmp(inputStreams{iI}, {remoteOutputs.name});
                
                if ~any(rI)
                    aas_log(aap, 1, sprintf('%s''s input stream ''%s'' does not come from any module in this AA, or from one of your remote locations.\nTry connecting the AA pipelines *after* all aas_addinitialstream() calls in your user script.', mod.name, inputStreams{iI}));
                end
                
                if sum(rI) > 1
                    aas_log(aap, 1, sprintf('%s is present at more than one remote location? Not sure how this happened', inputStreams{iI}));
                end
                
                remoteOutput = remoteOutputs(rI);
                remoteModule = remoteAA{remoteOutput.locI}.tasklist.main.module(remoteOutput.modI);
                remoteStreams(end+1) = struct('stream',       inputStreams{iI}, ...
                    'stagetag',     aas_getstagetag(remoteAA{remoteOutput.locI}, remoteOutput.modI), ...
                    'sourcedomain', remoteAA{remoteOutput.locI}.schema.tasksettings.(remoteModule.name)(remoteModule.index).ATTRIBUTE.domain, ...
                    'host',         remoteAAlocations(remoteOutput.locI).host, ...
                    'aapfilename',  fullfile(remoteAAlocations(remoteOutput.locI).directory, 'aap_parameters.mat'), ...
                    'allowcache',   remoteAAlocations(remoteOutput.locI).allowcache);
                
            end
            
        end
        
        aap.tasklist.main.module(modI).remotestream = remoteStreams;
        aap.aap_beforeuserchanges.tasklist.main.module(modI).remotestream = remoteStreams;
    end
    
    if isfield(aap.tasksettings.(mod.name)(mod.index), 'outputstreams')
        if isfield(aap.tasksettings.(mod.name)(mod.index).outputstreams, 'stream')
            outputStreams = aap.tasksettings.(mod.name)(mod.index).outputstreams.stream;
            if ~iscell(outputStreams), outputStreams = {outputStreams}; end
            
            % Update outputs present in the local AA
            for oI = 1 : length(outputStreams)
                if ~ismember(outputStreams{oI}, prevOutputs)
                    prevOutputs{end+1} = outputStreams{oI};
                end
            end
        end
    end
end


end

