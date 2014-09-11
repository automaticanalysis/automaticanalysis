function aap = aas_connectAApipelines(aap, remoteAAlocations)
%
% aas_connectAApipelines(aap, remoteAAlocations)
%
%
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
%
% -------------------------------------------------------
% Input parameters:
%
%   aap               : the AA structure
%
%   remoteAAlocations : a struct array with the following fields:
%
%     'host':        The IP address of the remote machine. Leave this as an
%                    empty string '' if the remote location is on the local
%                    machine.  SSH keys must be set up in advance to use
%                    this feature (see the WIKI).
%
%     'directory':   The full path to remote the AA pipeline. In that
%                    folder there should be file called aap_parameters.mat
%                    that contains all the AA information about the remote
%                    analysis.
%
%     'allowcache':  If set to 1, streams from the remote location will be
%                    cached on the local machine (in ~/aaworker/). Useful
%                    for preventing excessive network activity if you are
%                    connecting to remote machines.
%                    If set to 0, disables caching. Should probably be set
%                    to 0 if your remote AAPs live on the same local
%                    machine.
%                    If set to -1, will use caching behaviour determing by
%                    aap.directory_conventions.allowremotecache.
%
%     'maxstagetag': By default, we take data streams from the remote AAP
%                    from the *last* stage at which they are outputs.
%                    However, if you want to take stream outputs from an
%                    earlier stage, specify the full stage tag here.  Leave
%                    as an empty string '' for default behaviour, or give a
%                    string (e.g., 'aamod_realign_00001').
%
%     'checkMD5':    If set to 1, we will compare the status of all remote
%                    streams to the local versions before running AA
%                    pipeline. This way, we can detect if the remote data
%                    has changed in any way and cause local stages
%                    dependent on those data to reset. If set to 0, we
%                    don't check the remote data - once local modules that
%                    use those data are completed, they will never check
%                    again.
%
%     remoteAAlocations is struct array, so we can specify multiple
%     remote locations.  Priority is given to earlier array elements, so if
%     a data stream is found at multiple remote locations, we take that
%     stream from the location that occurs first in the array.
%
% -------------------------------------------------------
% created by:  cwild 2014-03-10
%
% updates:
%
% rhodri & cwild 2014-09-*: Update to allow fully qualified stream names in
% local and remote analyses. E.g., using aamod_realign_00001.epi to fetch
% the epi stream from the realign stage of a remote analysis, instead of
% the last occurence of epi.
% cwild 2014-09-09: output/input stream searching respects branches in the
% local and remote analyses.
% cwild 2014-04-02: Major update, added check for udpated data on the
% remote
% cwild 2014-03-18: misc cleaning
%

% Error checking:
if isempty(aap.acq_details.subjects)
    aas_log(aap, 1, 'aas_connectAApipelines() should be used after you have added subjects in your user script.');
end

if isempty(aap.acq_details.sessions)
    aas_log(aap, 1, 'aas_connectAApipelines() should be used after you have added sessions in your user script.');
end

if any(~isfield(remoteAAlocations, {'host', 'directory', 'allowcache', 'maxstagetag', 'checkMD5'}))
    aas_log(aap, 1, 'remoteAAlocations (input to aas_connectAApipelines) should be a struct array with the following fields: ''host'', ''directory'', ''allowcache'', ''maxstagetag'' ''checkMD5''');
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
remoteOutputs = struct('name', {{}}, 'locI', {}, 'modI', {}, 'stagetag', {});

local2remoteMaps = {};

% Collect output streams from the remote AA locations
% Reverse order, because priority is given to earlier array elements
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
    %
    %       % NOT BOTHERING TO CHECK SUBJECTS AND SESSIONS, BECAUSE WE
    %       DON'T WANT TO FORCE ALL DOMAINS TO BE PRESENT IN ALL AAs
    %
    %     % Check that all subjects in the local AA are present in the remote AA
    %     localSub = {aap.acq_details.subjects.mriname};
    %     remoteSub = {remoteAA{locI}.acq_details.subjects.mriname};
    %     subMatch = ismember(localSub, remoteSub);
    %     if any(~subMatch)
    %         aas_log(aap, 1, sprintf('Remote AAP (%s:%s) doesn''t have subjects: %s', remoteAAlocations(locI).host, remoteAAPfn, strjoin(localSub(~subMatch))));
    %     end
    %
    %     % Check that all sessions in the local AA are present in the remote AA
    %     if ~(length(aap.acq_details.sessions)==1 && isempty(aap.acq_details.sessions(1).name))
    %         localSess = {aap.acq_details.sessions.name};
    %         remoteSess = {remoteAA{locI}.acq_details.sessions.name};
    %         sessMatch = ismember(localSess, remoteSess);
    %         if any(~sessMatch)
    %             aas_log(aap, 1, sprintf('Remote AAP (%s:%s) doesn''t have sessions: %s', remoteAAlocations(locI).host, remoteAAPfn, strjoin(localSess(~sessMatch))));
    %         end
    %     end
    
    
    local2remoteMaps{locI} = aas_mapindices_betweenAAPs(aap, remoteAA{locI});
    
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
        
        if isfield(remoteAA{locI}.tasksettings.(mod.name)(mod.index), 'outputstreams') && ...
                isfield(remoteAA{locI}.tasksettings.(mod.name)(mod.index).outputstreams, 'stream')
            
            outputStreams = remoteAA{locI}.tasksettings.(mod.name)(mod.index).outputstreams.stream;

            if ~iscell(outputStreams), outputStreams = {outputStreams}; end
            
            % add all outputs to the end of the list
            for oI = 1 : length(outputStreams)
                remoteOutputs(end+1) = struct('name', outputStreams{oI}, 'locI', locI, 'modI', modI, 'stagetag', aas_getstagetag(remoteAA{locI}, modI));
            end
        end
    end
end

% Track the names of output streams from modules in the local analysis.  If
% a stream is an output from a previous stage in the local AA, then we
% don't bother trying to connect it from the remote AA.
allPrevOutputs = struct('streamname', {}, 'stagetag', {}, 'moduleInd', {}, 'dependentOn', {});

aas_log(aap,0,'Checking status of remote streams...');

for modI = 1 : length(aap.tasklist.main.module)
    mod = aap.tasklist.main.module(modI);
    stagetag = aas_getstagetag(aap, modI);
    
    if isfield(aap.tasksettings.(mod.name)(mod.index),'inputstreams') && ...
            isfield(aap.tasksettings.(mod.name)(mod.index).inputstreams, 'stream')
        
        % Names of input and output streams for this module
        inputStreams = aap.tasksettings.(mod.name)(mod.index).inputstreams.stream;
        if ~iscell(inputStreams), inputStreams = {inputStreams}; end
        
        remoteStreams = struct('stream', {}, 'stagetag', {}, 'sourcedomain', {}, 'host', {}, 'aapfilename', {}, 'allowcache', {});
        
        for iI = 1 : length(inputStreams)

            % Is it fully qualified, with source module stage tag as well (e.g.,
            %   <stream>aamod_coreg_extended_1_00001.structural</stream>
            if any(inputStreams{iI} == '.')
                [inputStageTag, rem] = strtok(inputStreams{iI}, '.');
                inputStreamName = strtok(rem,'.');
            else
                inputStageTag = '';
                inputStreamName = inputStreams{iI};
            end
            
            % If the input doesn't come from a previous module, let's add it to
            % the list of remote streams for this module.
            if isempty(allPrevOutputs)
                foundPrevious=false;
            else
                
                if ~isempty(inputStageTag)
                    foundPrevious = any([strcmp(inputStreamName,{allPrevOutputs.streamname})] & [strcmp(inputStageTag,{allPrevOutputs.stagetag})]);
                else
                    
                    % not fully qualified, follow the curent branch back
                    % up to the top looking for the stream
                    foundPrevious = false;
                    prevModI = -1;
                    if ~isempty(mod.tobecompletedfirst)
                        prevModI = aas_getmoduleindexfromtag(aap, mod.tobecompletedfirst{1});
                    end
                    while ~foundPrevious && prevModI > 0
                        prevModOutputsI = [allPrevOutputs.moduleInd] == prevModI; % outputs of the module that comes before this one in this branch
                        prevModOutputs = allPrevOutputs(prevModOutputsI);
                        
                        if isempty(prevModOutputs)
                            prevModI = aas_getmoduleindexfromtag(aap, aap.tasklist.main.module(prevModI).tobecompletedfirst{1});
                        else
                            foundPrevious = ismember(inputStreamName, {prevModOutputs.streamname});
                            prevModI = prevModOutputs(1).dependentOn;
                        end
                    end
                end
            end
            
            if ~foundPrevious
                % Check previously added remotestreams
                if ~isempty(mod.remotestream) && ismember(inputStreams{iI}, {mod.remotestream.stream}), continue; end
                
                % Is the input stream present in the list of remote streams?
                rI = find(strcmp(inputStreamName, {remoteOutputs.name}));
                
                % Empty stagetag means non-qualified name, take last
                % occurence of this stream
                if ~isempty(inputStageTag) && any(rI)

                    rIQual = find(strcmp(inputStageTag, {remoteOutputs(rI).stagetag}));
                    
                    if isempty(rIQual)
                        aas_log(aap, 1, sprintf('Can''t find %s.%s in remotely specified AAPs!', inputStageTag, inputStreamName));
                    end
                    rI = rI(rIQual(end));
                end
                
                if any(rI)
                    
                    rI = rI(end);
                    
                    remoteOutput = remoteOutputs(rI);
                    remoteModule = remoteAA{remoteOutput.locI}.tasklist.main.module(remoteOutput.modI);
                    
                    trgDomain = aap.schema.tasksettings.(mod.name)(mod.index).ATTRIBUTE.domain;
                    srcDomain = remoteAA{remoteOutput.locI}.schema.tasksettings.(remoteModule.name)(remoteModule.index).ATTRIBUTE.domain;
                    
                    remoteStreams(end+1) = struct('stream',       inputStreams{iI}, ...
                        'stagetag',     aas_getstagetag(remoteAA{remoteOutput.locI}, remoteOutput.modI), ...
                        'sourcedomain', srcDomain, ...
                        'host',         remoteAAlocations(remoteOutput.locI).host, ...
                        'aapfilename',  fullfile(remoteAAlocations(remoteOutput.locI).directory, 'aap_parameters.mat'), ...
                        'allowcache',   remoteAAlocations(remoteOutput.locI).allowcache);
                    
                    % Now check if the remote data has changed...
                    if remoteAAlocations(remoteOutput.locI).checkMD5
                        
                        % First we have to find all possible locations of the target streams
                        trgDomainTree = aas_dependencytree_finddomain(trgDomain, aap.directory_conventions.parallel_dependencies, {});
                        NPerDomain = [];
                        for tdI = 1 : length(trgDomainTree)
                            NPerDomain(tdI) = aas_getN_bydomain(aap, trgDomainTree{tdI}, []);
                        end
                        
                        % Generate all possible indices for the target domain
                        trgIndices = zeros(prod(NPerDomain), length(trgDomainTree));
                        for tdI = 1: length(trgDomainTree)
                            trgIndices(:,tdI) = repmat(kron([1:NPerDomain(tdI)]', ones(prod(NPerDomain(tdI+1:end)), 1)), prod(NPerDomain(1:tdI-1)), 1);
                        end
                        
                        trgDomainTree = trgDomainTree(2:end);
                        trgIndices(:,1) = [];
                        
                        % Check doneflags / MD5s for each occurence of this stream.
                        for tdI = 1 : size(trgIndices, 1)
                            
                            resetThisStage = false;
                            
                            % If this stage has completed, but target data has changed, let's force it to re-run
                            trgDoneFlag = aas_doneflag_getpath_bydomain(aap, trgDomain, trgIndices(tdI,1:find(strcmp(trgDomain, trgDomainTree))), modI);
                            if (aas_doneflagexists(aap, trgDoneFlag))
                                
                                % It's also possible for streams to be written
                                % to any subdomain of the source domain so we
                                % have to go hunt for all of those.  This makes
                                % this complicated :(
                                possibleLocs = aas_getdependencies_bydomain(aap, srcDomain, trgDomain, trgIndices(tdI,1:find(strcmp(trgDomain, trgDomainTree))));
                                
                                for lI = 1 : length(possibleLocs)
                                    
                                    % Stream header file with MD5
                                    trgPath = aas_getpath_bydomain(aap, possibleLocs{lI}{1}, possibleLocs{lI}{2}, modI);
                                    trgStreamDesc = fullfile(trgPath, sprintf('stream_%s_inputto_%s.txt', inputStreams{iI}, sprintf('%s_%05d', mod.name, mod.index)));
                                    
                                    if exist(trgStreamDesc)
                                        
                                        % Map the local indices to the remote AA
                                        localIndices = possibleLocs{lI}{2};
                                        remoteIndices = [];
                                        srcDomainTree = aas_dependencytree_finddomain(possibleLocs{lI}{1}, remoteAA{locI}.directory_conventions.parallel_dependencies, {});
                                        srcDomainTree = srcDomainTree(2:end);
                                        for i = 1 : length(localIndices)
                                            try remoteIndices(i) = local2remoteMaps{remoteOutput.locI}.(srcDomainTree{i})(localIndices(i));
                                            catch
                                                disp('');
                                            end
                                        end
                                        
                                        % A remote index of 0 means that we can't
                                        % find that one, so it just might not be
                                        % present in the remote location.  That's
                                        % OK because maybe we're pulling it from
                                        % somewhere else....
                                        if ~any(remoteIndices==0)
                                            
                                            % Copy over the stream desc from the remote
                                            remoteStage = aas_getstagetag(remoteAA{locI},remoteOutput.modI);
                                            try remoteSrcPath = aas_getpath_bydomain(remoteAA{remoteOutput.locI}, possibleLocs{lI}{1}, remoteIndices, remoteOutput.modI);
                                            catch
                                                disp('');
                                            end
                                            remoteSrcDesc = fullfile(remoteSrcPath, sprintf('stream_%s_outputfrom_%s.txt', strrep(inputStreams{iI},[remoteStage '.'],''), remoteStage));
                                            localSrcDesc = fullfile(trgPath, sprintf('stream_%s_remoteoutputfrom_%s_%s.txt', inputStreams{iI}, remoteAAlocations(remoteOutput.locI).host, remoteStage));
                                            aap = aas_copyfromremote(aap, remoteAAlocations(remoteOutput.locI).host, remoteSrcDesc, localSrcDesc, 'allow404', 1, 'allowcache', remoteAAlocations(remoteOutput.locI).allowcache, 'verbose', 0);
                                            
                                            % If it didn't copy, then we should reset this stage
                                            if ~exist(localSrcDesc)
                                                resetThisStage = true;
                                                aas_log(aap, 0, sprintf('%s %s input %s REMOTE status UNKNOWN. Forcing this stage to re-run.', mod.name, strjoin(arrayfun(@(x) sprintf('[%d]',x), trgIndices(tdI,:), 'UniformOutput', false)), remoteOutput.name), 'red');
                                            else
                                                
                                                % Compare the MD5s
                                                trgFileID = fopen(trgStreamDesc, 'r');
                                                [aap, trgMD5] = aas_load_md5(aap, trgFileID, remoteOutput.name);
                                                
                                                srcFileID = fopen(localSrcDesc, 'r');
                                                [aap, srcMD5] = aas_load_md5(aap, srcFileID, remoteOutput.name);
                                                
                                                fclose(trgFileID); fclose(srcFileID);
                                                
                                                if ~strcmp(trgMD5, srcMD5)
                                                    resetThisStage = true;
                                                    aas_log(aap, 0, sprintf('%s %s input %s REMOTE status CHANGED. Forcing this stage to re-run.', mod.name, strjoin(arrayfun(@(x) sprintf('[%d]',x), trgIndices(tdI,:), 'UniformOutput', false)), remoteOutput.name), 'red');
                                                end
                                            end
                                        end
                                    end
                                    
                                    if resetThisStage
                                        aas_delete_doneflag_bypath(aap, trgDoneFlag);
                                    end
                                    
                                end % End if done flag exists
                                
                            end % End check target destination indices
                            
                        end % End for each possible stream location
                        
                    end % End if check MD5s
                    
                else
                    
                    if isstruct(aap.schema.tasksettings.(mod.name)(mod.index).inputstreams.stream{iI}) && ...
                            isfield(aap.schema.tasksettings.(mod.name)(mod.index).inputstreams.stream{iI}, 'isessential') && ...
                            aap.schema.tasksettings.(mod.name)(mod.index).inputstreams.stream{iI}.ATTRIBUTE.isessential
                        aas_log(aap, 1, sprintf('%s''s input stream ''%s'' does not come from any module in this AA, or from one of your remote locations.\nTry connecting the AA pipelines *after* all aas_addinitialstream() calls in your user script.', mod.name, inputStreams{iI}));
                    end
                    
                end
            end % End if input ~present in prev outputs
            
        end % End loop over input straems
        
        aap.tasklist.main.module(modI).remotestream = remoteStreams;
        aap.aap_beforeuserchanges.tasklist.main.module(modI).remotestream = remoteStreams;
    end
    
    if isfield(aap.tasksettings.(mod.name)(mod.index), 'outputstreams')
        if isfield(aap.tasksettings.(mod.name)(mod.index).outputstreams, 'stream')
            outputStreams = aap.tasksettings.(mod.name)(mod.index).outputstreams.stream;
            if ~iscell(outputStreams), outputStreams = {outputStreams}; end
            
            % Update outputs present in the local AA
            for oI = 1 : length(outputStreams)
                allPrevOutputs(end+1) = struct('streamname', outputStreams{oI}, 'stagetag', stagetag, 'moduleInd', modI, 'dependentOn', -1);
                if ~isempty(mod.tobecompletedfirst)
                    allPrevOutputs(end).dependentOn = aas_getmoduleindexfromtag(aap, mod.tobecompletedfirst{1});
                end
            end
        end
    end
end


end


