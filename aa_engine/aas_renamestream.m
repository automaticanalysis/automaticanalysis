function aap=aas_renamestream(aap,stage,origstream,newstream,type)
% AAS_RENAMESTREAM 
% Rename a TYPE (input|output) of stream of a STAGE from ORIGSTREAM to NEWSTREAM. 
% If ORIGSTREAM is 'append', then the NEWSTREAM will be appended to the list of streams
% If NEWSTREAM is [], then the ORIGSTREAM will be removed from the list of streams
%
% E.g. aap = aas_renamestream(aap,'aamod_firstlevel_model_00001','epi','aamod_coreg_extended_2epi_00001.epi')

if nargin < 5
    type = 'input';
end

if strcmp(type,'output')
    aas_log(aap,false,sprintf('WARNING: Renaming output %s of %s to %s.\n  It might breake dependency for the consecuting module(s).\n  Make sure that you know what you are doing!',origstream,stage,newstream));
end

%% Locate STAGE
[stagename, index] = strtok_ptrn(stage,'_0');
stageindex = sscanf(index,'_%05d');
if ~isfield(aap.tasksettings,stagename) || (numel(aap.tasksettings.(stagename)) < stageindex)
    aas_log(aap,1,sprintf('ERROR: Stage %s not found!',stage));
end
modulenumber = cell_index({aap.tasklist.main.module.name},stagename); modulenumber = modulenumber(stageindex);

%% Check task
odot = origstream == '.';
if any(odot), origstreamname = origstream(find(odot)+1:end);
else, origstreamname = origstream; end

if ~isempty(newstream)
    dat = strsplit(newstream,':');
    newstream = dat{1}; newattr = dat(2:end);
else
    newattr = '';
end
ndot = newstream == '.';
if any(ndot), newstreamname = newstream(find(ndot)+1:end);
else, newstreamname = newstream; end

toSpecify = any(ndot);
toRename = ~strcmp(origstreamname,newstreamname);

%% Locate ORIGSTREAM
inputstreams = aap.tasksettings.(stagename)(stageindex).([type 'streams']).stream;
if ~iscell(inputstreams)
    inputstreams = {inputstreams};
end
for i = 1:numel(inputstreams)
    idot = inputstreams{i} == '.';
    if any(idot), inputstreamname = inputstreams{i}(find(idot)+1:end);
    else, inputstreamname = inputstreams{i}; end
    if strcmp(inputstreamname,origstreamname)
        break;
    end
end

stream = aap.schema.tasksettings.(stagename)(stageindex).([type 'streams']).stream{i};
if strcmp(inputstreamname,origstreamname)
    ind = i;
elseif strcmp(origstreamname,'append')
    ind = i + 1;    
    stream = [];
    stream.CONTENT = newstream;
    % defaults
    stream.ATTRIBUTE.isrenameable = 1;
    stream.ATTRIBUTE.isessential = 1;
    for attr = newattr
        keyval = strsplit(attr{1},'-');
        stream.ATTRIBUTE.(keyval{1}) = str2double(keyval{2});
    end
    aap.schema.tasksettings.(stagename)(stageindex).([type 'streams']).stream{ind} = stream;
    if isfield(aap,'internal'), aap.internal.aap_initial.schema.tasksettings.(stagename)(stageindex).([type 'streams']).stream{ind} = stream; end
else
    aas_log(aap,1,sprintf('ERROR: Stream %s of module %s not found!',origstream,stagename));
end

switch type
    case 'input'
        internalstore = 'inputstreamsources';
    case 'output'
        internalstore = 'outputstreamdestinations';
end
if isfield(aap,'internal')
    storedstreams = aap.internal.(internalstore){modulenumber}.stream;
    if isempty(storedstreams) || ~any(cell_index({storedstreams.name},origstreamname)), streamindex = NaN;
    else
        if ~strcmp(origstreamname,'append'), streamindex = cell_index({storedstreams.name},origstreamname); 
        else, streamindex = numel(storedstreams); end
    end
end
        
%% Check ORIGSTREAM
if (~isstruct(stream) || ~isfield(stream.ATTRIBUTE,'isrenameable') || ~stream.ATTRIBUTE.isrenameable) && toRename % no attributes --> assume non-renameable --> allow specifying only
    aas_log(aap,1,sprintf('ERROR: Stream %s of module %s does not support renaming!\nERROR: Check attribute "isrenameable"!',origstream,stagename));
end
if toSpecify
    sourcestage = newstream(1:find(ndot)-1);
    [sourcestagename, index] = strtok_ptrn(sourcestage,'_0');
    sourcestageindex = sscanf(index,'_%05d');
    if ~isfield(aap.tasksettings,sourcestagename) || (numel(aap.tasksettings.(sourcestagename)) < sourcestageindex)
        aas_log(aap,0,sprintf('WARNING: Stage %s not found in local pipeline!',sourcestage));
    end
end

%% Rename ORIGSTREAM
if ~iscell(aap.tasksettings.(stagename)(stageindex).([type 'streams']).stream) % single stream
    singlenewstream = newstream;
    if strcmp(origstreamname,'append'), singlenewstream = {aap.tasksettings.(stagename)(stageindex).([type 'streams']).stream newstream}; end
    aap.tasksettings.(stagename)(stageindex).([type 'streams']).stream = singlenewstream;
    aap.aap_beforeuserchanges.tasksettings.(stagename)(stageindex).([type 'streams']).stream = singlenewstream;
    if isfield(aap,'internal')
        aap.internal.aap_initial.tasksettings.(stagename)(stageindex).([type 'streams']).stream = singlenewstream;
        aap.internal.aap_initial.aap_beforeuserchanges.tasksettings.(stagename)(stageindex).([type 'streams']).stream = singlenewstream;
    end
else
    aap.tasksettings.(stagename)(stageindex).([type 'streams']).stream{ind} = newstream;
    aap.aap_beforeuserchanges.tasksettings.(stagename)(stageindex).([type 'streams']).stream{ind} = newstream;
    if isfield(aap,'internal')
        aap.internal.aap_initial.tasksettings.(stagename)(stageindex).([type 'streams']).stream{ind} = newstream;
        aap.internal.aap_initial.aap_beforeuserchanges.tasksettings.(stagename)(stageindex).([type 'streams']).stream{ind} = newstream;
    end
end

% Edit internalstore (inputstreamsources or outputstreamdestinations)
if isfield(aap,'internal') && ((numel(streamindex) > 1) || ~isnan(streamindex))
    if isempty(newstream) % remove
        aap.internal.(internalstore){modulenumber}.stream(streamindex) = [];
        aap.internal.aap_initial.internal.(internalstore){modulenumber}.stream(streamindex) = [];
    else
        for i = streamindex'
            aap.internal.(internalstore){modulenumber}.stream(i).name = newstream;
            aap.internal.aap_initial.internal.(internalstore){modulenumber}.stream(i).name = newstream;
        end
    end
end

% Edit schema
if isempty(newstream) || ~isstruct(aap.schema.tasksettings.(stagename)(stageindex).([type 'streams']).stream{ind})
    aap.schema.tasksettings.(stagename)(stageindex).([type 'streams']).stream{ind} = newstream;
    if isfield(aap,'internal'), aap.internal.aap_initial.schema.tasksettings.(stagename)(stageindex).([type 'streams']).stream{ind} = newstream; end
elseif ~isempty(newstream) % with attributes
    aap.schema.tasksettings.(stagename)(stageindex).([type 'streams']).stream{ind}.CONTENT = newstream;
    if isfield(aap,'internal'), aap.internal.aap_initial.schema.tasksettings.(stagename)(stageindex).([type 'streams']).stream{ind}.CONTENT = newstream; end
end

% Update current (if applicable)
if isfield(aap.tasklist,'currenttask'), aap.tasklist.currenttask.([type 'streams']) = aap.schema.tasksettings.(stagename)(stageindex).([type 'streams']); end
