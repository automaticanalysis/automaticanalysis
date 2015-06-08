function aap=aas_renamestream(aap,stage,origstream,newstream,type)
% AAS_RENAMESTREAM 
% Rename a TYPE (input|output) of stream of a STAGE from ORIGSTREAM to NEWSTREAM

if nargin < 5
    type = 'input';
end

if strcmp(type,'output')
    aas_log(aap,true,'Modifying outputstreams is not implemented yet!');
end

%% Locate STAGE
[stagename, index] = strtok_ptrn(stage,'_0');
index = sscanf(index,'_%05d');
stages = {aap.tasklist.main.module.name};
ind = find(strcmp(stages,stagename),index,'first');
if isempty(ind)
    aas_log(aap,1,sprinf('Stage %s not found!',stage));
end
stageindex = ind(end);

%% Locate ORIGSTREAM
ind = cell_index(aap.tasksettings.(stagename)(end).([type 'streams']).stream,origstream);

%% Rename ORIGSTREAM
aap.tasksettings.(stagename)(end).([type 'streams']).stream{ind}=newstream;
aap.aap_beforeuserchanges.tasksettings.(stagename)(end).([type 'streams']).stream{ind}=newstream;
aap.schema.tasksettings.(stagename)(end).([type 'streams']).stream{ind}=newstream;
