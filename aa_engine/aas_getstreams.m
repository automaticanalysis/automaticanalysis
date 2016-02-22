function [streams, attributes] = aas_getstreams(varargin)
aap = varargin{1};
select = varargin{end};

if nargin == 2 % current
    task = aap.tasklist.currenttask;
    [stagename, ind] = strtok_ptrn(task.name,'_0');
    index = sscanf(ind,'_%d');
else % any
    stagename = varargin{2};
    index = 1;
    if nargin == 4, index = varargin{3}; end
    task = aap.tasksettings.(stagename)(index);
end

streams = task.([select 'streams']).stream;
if ~iscell(streams), streams = {streams}; end
schemas = aap.schema.tasksettings.(stagename)(index).([select 'streams']).stream;

% Remove/save attributes
attributes = {};
for s = 1:numel(streams)
    if isstruct(streams{s}), streams{s} = streams{s}.CONTENT; end
    if isstruct(schemas{s}), attributes{s} = schemas{s}.ATTRIBUTE; end % assume same order
end