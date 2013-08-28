%%==============================
% Delete done flags
% function aas_delete_doneflag(aap,stage,domain,indices)

function aap=aas_delete_doneflag_bydomain(aap,stage,outputdomain,indices)

% allow full path of module to be provided
[stagepath stagename]=fileparts(aap.tasklist.main.module(stage).name);
index=aap.tasklist.main.module(stage).index;
% find out whether this module needs to be executed once per study, subject or session
domain=aap.schema.tasksettings.(stagename)(index).ATTRIBUTE.domain;

% The following was written to find where the inputs to a stage come from.
% Here, we want to know where the outputs of a stage are going. So, the
% "target" and "source" domains are reversed
deps=aas_getdependencies_bydomain(aap,domain,outputdomain,indices,'doneflaglocations_thatexist',stage);
for depind=1:length(deps)
    doneflag=aas_doneflag_getpath_bydomain(aap,deps{depind}{1},deps{depind}{2},stage);
    aas_delete_doneflag_bypath(aap,doneflag);
end;
