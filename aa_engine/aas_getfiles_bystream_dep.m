function fname = aas_getfiles_bystream_dep(aap,domain,indices,stream)
% locate inputstream
inputstreamindex = strcmp({aap.internal.inputstreamsources{aap.tasklist.currenttask.modulenumber}.stream.name},stream);
if not(any(inputstreamindex))
    inputstreamindex = cell2mat(cellfun(@(x) any(strcmp(strsplit(x,'.'),stream)),{aap.internal.inputstreamsources{aap.tasklist.currenttask.modulenumber}.stream.name},'UniformOutput',false));
end

% obnain source module number
sourcenumber = aap.internal.inputstreamsources{aap.tasklist.currenttask.modulenumber}.stream(inputstreamindex).sourcenumber;

% switch source module
aap = aas_setcurrenttask(aap,sourcenumber);

% obtain file
fname = aas_getfiles_bystream_multilevel(aap,domain,indices,stream,'output');
end