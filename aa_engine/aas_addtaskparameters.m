function [aap index mfile_alias]=aas_addtaskparameters(aap,taskname_full,tasknamealias_full)

% Can provide an alias, which is used for the .xml name
if (~exist('tasknamealias_full','var') || isempty(tasknamealias_full))
    tasknamealias_full=taskname_full;
end;

[taskpath taskname]=fileparts(taskname_full); % allow full path to be passed
[taskpathalias tasknamealias]=fileparts(tasknamealias_full); % allow full path to be passed

addpath(taskpathalias,'-END'); % if full path was provided, xml file probably in same location
% module can also now be called by name

xmlmodfn=which([tasknamealias '.xml']);
if isempty(xmlmodfn)
    aas_log(aap,true,sprintf('Each module needs an accompanying XML header, could not find one for %s.',taskname));
end;
Pref.ReadAttr=0;
xml_module=xml_read(xmlmodfn,Pref);
xml_module_schema=xml_read(xmlmodfn);

try
    tasksettings=xml_module.tasksettings.(taskname);
    tasksettings_schema=xml_module_schema.tasksettings.(taskname);
catch
    try
        tasksettings=xml_module.tasklist.currenttask;
        tasksettings_schema=xml_module_schema.tasklist.currenttask;
    catch
        aas_log(aap,true,sprintf('Was expecting XML file %s to have branch either \naap.tasklist.currenttask\nor\naap.tasksettings.%s',xmlmodfn,taskname));
    end;
end;

% Added ability to inherit from other .xml module headers [CW 2014-01-14]
[tasksettings tasksettings_schema] = recursiveInheritFrom(tasksettings, tasksettings_schema);

% Put in mfile_alias if we're using an aliased name, and no mfile_alias
% already provideid
mfile_alias='';
if (~strcmp(taskname_full,tasknamealias_full))
    if (~isfield(tasksettings_schema.ATTRIBUTE,'mfile_alias'))
        tasksettings_schema.ATTRIBUTE.mfile_alias=tasknamealias_full;
    end;
end;

% if this is the first task to be specified, add tasksettings field to aap
if ~isfield(aap,'tasksettings'); aap.tasksettings=struct([]); end

if (isfield(aap.tasksettings,taskname))
    index=length(aap.tasksettings.(taskname))+1;
    %    fprintf('Task %s index %d\n',taskname,index);
    if ~isfield(tasksettings_schema.ATTRIBUTE,'startempty')
        tmp=tasksettings;
        tmp.timeadded=clock;
        aap.tasksettings.(taskname)(index)=tmp;
    else
        aap.tasksettings.(taskname)(index)= struct('desc',tasksettings_schema.ATTRIBUTE.desc);
    end
else
    if ~isfield(tasksettings_schema.ATTRIBUTE,'startempty')
        tmp=tasksettings;
        tmp.timeadded=clock;
        aap.tasksettings(1).(taskname)=tmp;
    else
        aap.tasksettings(1).(taskname).desc=tasksettings_schema.ATTRIBUTE.desc;
    end
    index=1;
    %    fprintf('Task %s index %d [start]\n',taskname,index);
end;



% ...and update before user changes to reflect new task options
aap.aap_beforeuserchanges.tasksettings.(taskname)(index)=aap.tasksettings.(taskname)(index);

% And update schema

% Allow for empty inputstreams
if ~isfield(tasksettings_schema,'inputstreams')
    tasksettings_schema.inputstreams=[];
end;
if ~isfield(tasksettings_schema.inputstreams,'stream')
    tasksettings_schema.inputstreams.stream={};
end;
% Turn single streams into cells
if ~iscell(tasksettings_schema.inputstreams.stream)
    tasksettings_schema.inputstreams.stream={tasksettings_schema.inputstreams.stream};
end;
% Turn single cell of array of structs into cell of structs
if length(tasksettings_schema.inputstreams.stream)==1 && ...
        isstruct(tasksettings_schema.inputstreams.stream{1}) && ...
        isfield(tasksettings_schema.inputstreams.stream{1},'ATTRIBUTE')
    streams = tasksettings_schema.inputstreams.stream{1};
    for i = 1:numel(streams)
        tasksettings_schema.inputstreams.stream{i} = streams(i);
    end
end;

% Allow for empty outputstreams
if ~isfield(tasksettings_schema,'outputstreams')
    tasksettings_schema.outputstreams=[];
end;
if ~isfield(tasksettings_schema.outputstreams,'stream')
    tasksettings_schema.outputstreams.stream={};
end;
% Turn single streams into cells
if ~iscell(tasksettings_schema.outputstreams.stream)
    tasksettings_schema.outputstreams.stream={tasksettings_schema.outputstreams.stream};
end;
% Turn single cell of array of structs into cell of structs
if length(tasksettings_schema.outputstreams.stream)==1 && ...
        isstruct(tasksettings_schema.outputstreams.stream{1}) && ...
        isfield(tasksettings_schema.outputstreams.stream{1},'ATTRIBUTE')
    streams = tasksettings_schema.outputstreams.stream{1};
    for i = 1:numel(streams)
        tasksettings_schema.outputstreams.stream{i} = streams(i);
    end
end;

aap.schema.tasksettings.(taskname)(index)=tasksettings_schema;

% ...and update before user changes to reflect new task options
aap.aap_beforeuserchanges.schema.tasksettings.(taskname)(index)=tasksettings_schema;

end

function [tasksettings, tasksettings_schema] = recursiveInheritFrom(tasksettings, tasksettings_schema)

% If we inherit settings from another module....
if isfield(tasksettings_schema.ATTRIBUTE, 'inheritfrom') && ~isempty(tasksettings_schema.ATTRIBUTE.inheritfrom)
    
    % Find the parent module .xml
    inheritFrom = [tasksettings_schema.ATTRIBUTE.inheritfrom '.xml'];
    if ~exist(inheritFrom), error('Can''t find %s', inheritFrom); end
    
    % Get the parent's settings
    Pref.ReadAttr=0;
    parentSettings = xml_read(inheritFrom, Pref);
    parentSettings = parentSettings.tasklist.currenttask;
    
    % get the parent's schema
    parentSettings_schema = xml_read(inheritFrom);
    parentSettings_schema = parentSettings_schema.tasklist.currenttask;
    
    % Does the parent inherit settings, too?
    [parentSettings parentSettings_schema] = recursiveInheritFrom(parentSettings, parentSettings_schema);
    
    % Now replace the settings in the parent module with those from
    % the child module, or add the setting if it wasn't there before.
    fields = fieldnames(tasksettings);
    for f = 1:numel(fields)
        parentSettings.(fields{f}) = tasksettings.(fields{f});
    end
    tasksettings = parentSettings;
    
    fields = fieldnames(tasksettings_schema);
    for f = 1:numel(fields)
        parentSettings_schema.(fields{f}) = tasksettings_schema.(fields{f});
    end
    tasksettings_schema = parentSettings_schema;
    
end

end