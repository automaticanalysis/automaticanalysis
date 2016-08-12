% Retrieve streamnames of a stage
%
% FORMAT [streams, attributes] = aas_getstreams(aap,[modulename,][moduleindex,]streamtype[,selection])
%   Outputs
%   - stream:       streamnames (cell in case of multiple streams)
%   - attributes:   stream attributes (cell in case of multiple streams)
%   Inputs
%   - aap:          aap structure of the current stage
%   - modulename:   modulename of the stage
%   - moduleindex:  moduleindex of the stage
%   - streamtype:   'input'|'output'
%   - selection:    index of the stream (can be a single number of an array)

function [streams, attributes] = aas_getstreams(aap,varargin)

%% Parse input
if any(strcmp(varargin{1},{'input','output'})) % current
    task = aap.tasklist.currenttask;
    [stagename, ind] = strtok_ptrn(task.name,'_0');
    index = sscanf(ind,'_%d');
else % any
    stagename = varargin{1}; varargin(1) = [];
    index = 1; if isnumeric(varargin{1}), index = varargin{1}; varargin(1) = []; end
    task = aap.tasksettings.(stagename)(index);
end
streamtype = varargin{1}; varargin(1) = [];

select = []; if ~isempty(varargin), select = varargin{1}; end

%% Process
streams = task.([streamtype 'streams']).stream;
if ~iscell(streams), streams = {streams}; end
schemas = aap.schema.tasksettings.(stagename)(index).([streamtype 'streams']).stream;

% Remove/save attributes
attributes = cell(size(streams));
for s = 1:1:numel(streams)
    if isstruct(streams{s}), streams{s} = streams{s}.CONTENT; end
    if isstruct(schemas{s}), attributes{s} = schemas{s}.ATTRIBUTE; end % assume same order
end

switch numel(select)
    case 0
    case 1
        streams = streams{select};
        attributes = attributes{select};
    otherwise
        streams = streams(select);
        attributes = attributes(select);
end