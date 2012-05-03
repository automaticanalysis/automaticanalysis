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

aap.schema.tasksettings.(taskname)(index)=tasksettings_schema;

% ...and update before user changes to reflect new task options
aap.aap_beforeuserchanges.schema.tasksettings.(taskname)(index)=tasksettings_schema;