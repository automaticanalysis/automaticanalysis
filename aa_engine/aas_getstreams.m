function streams = aas_getstreams(aap,select)
streams = aap.tasklist.currenttask.([select 'putstreams']).stream;
for s = 1:numel(streams)
    if isstruct(streams{s}), streams{s} = streams{s}.CONTENT; end
end