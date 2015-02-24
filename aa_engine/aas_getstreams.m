function streams = aas_getstreams(varargin)
aap = varargin{1};
select = varargin{end};

if nargin == 2 % current
    task = aap.tasklist.currenttask;
else % any
    stagename = varargin{2};
    index = 1;
    if nargin == 4, index = varargin{3}; end
    task = aap.tasksettings.(stagename)(index);
end

streams = task.([select 'putstreams']).stream;
if ~iscell(streams), streams = {streams}; end

% Remove attributes
for s = 1:numel(streams)
    if isstruct(streams{s}), streams{s} = streams{s}.CONTENT; end
end