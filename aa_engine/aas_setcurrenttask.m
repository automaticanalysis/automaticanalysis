% This is used to set up particular components of the aap structure that
% change from module to module. It implements the ability to provide
% module-specific parameters (e.g., for branched pipelines) by applying the
% values in extraparameters.aap to the aap structure after setting it up

function [aap]=aas_setcurrenttask(aap,k)

% Start with the initial pure aap
initaap=aap.internal.aap_initial;
aap=initaap;
aap.internal.aap_initial=initaap;

stagename=aas_getstagetag(aap,k);

% Stuff the task likes to know
aap.tasklist.currenttask.epiprefix=aap.tasklist.main.module(k).epiprefix;
aap.tasklist.currenttask.extraparameters=aap.tasklist.main.module(k).extraparameters;
aap.tasklist.currenttask.settings=aap.tasksettings.(aap.tasklist.main.module(k).name)(aap.tasklist.main.module(k).index);
aap.tasklist.currenttask.inputstreams=aap.schema.tasksettings.(aap.tasklist.main.module(k).name)(aap.tasklist.main.module(k).index).inputstreams;
aap.tasklist.currenttask.outputstreams=aap.schema.tasksettings.(aap.tasklist.main.module(k).name)(aap.tasklist.main.module(k).index).outputstreams;
aap.tasklist.currenttask.name=stagename;
aap.tasklist.currenttask.index=aap.tasklist.main.module(k).index;
aap.tasklist.currenttask.modulenumber=k;
aap.tasklist.currenttask.domain=aap.schema.tasksettings.(aap.tasklist.main.module(k).name)(aap.tasklist.main.module(k).index).ATTRIBUTE.domain;

% Recursively copy and parameters in extraparameters.aap for this task on
% top of aap if provided
if (isfield(aap.tasklist.main.module(k).extraparameters,'aap'))
    aap=aas_copyparameters(aap.tasklist.main.module(k).extraparameters.aap,aap,'aap');
end;


% Check the apparent study root is set appropriately
aap.acq_details.root=aas_getstudypath(aap,k);
% ..and for remote filesystem if we're using one
remotefilesystem=aap.directory_conventions.remotefilesystem;
if (~strcmp(remotefilesystem,'none'))
    aap.acq_details.(remotefilesystem).root=aas_getstudypath(aap,remotefilesystem,k);
else
    aap.acq_details.root=aas_getstudypath(aap,k);
end;

% Check subselected sessions
selected_sessions=aap.acq_details.selected_sessions;
if (ischar(selected_sessions))
    if (strcmp(selected_sessions,'*'))
        % Wildcard, same as empty
        selected_sessions=1:length(aap.acq_details.sessions);
    else
        % Named sessions, parse to get numbers
        rem=strtrim(selected_sessions);
        selected_sessions=[];
        while(length(rem)>0)
            [sessionnme rem]=strtok(rem,' ');
            sessionind=find(strcmp({aap.acq_details.sessions.name},sessionnme));
            if (isempty(sessionind))
                aas_log(aap,true,sprintf('Unknown session %s specified in selected_sessions field of a branch in the tasklist, sessions were %s',sessionnme,sprintf('%s ',aap.acq_details.sessions.name)));
            end;
            selected_sessions=[selected_sessions sessionind];
        end;
    end;
    aap.acq_details.selected_sessions=selected_sessions;
end;

    % Set SPM defaults appropriately
    global defaults
    defaults=aap.spm.defaults;
    
end


