% Automatic analysis - add a task to tasklist
% Called mainly by aap_tasklist but also by user scripts
% Third parameter is prefix for input EPIs if they are taken
% function [aap]=aas_addtask(aap,taskname[,epiprefix[,tobecompletedfirst]][,extraparameters])
% 
%  extraparameters (which is a structure) may be specified without other
%  preceding string parameters epiprefix and tobecompletedfirst
%
%  Examples of usage:    
%   aap=aas_addtask(aap,'aamod_smooth');
%
%   extraparameters=[];
%   extraparameters.aap.tasksettings.aamod_smooth=14;
%   aap=aas_addtask(aap,'aamod_smooth',extraparameters);
%
%   aap=aas_addtask(aap,'aamod_smooth','wuf','aamod_realign');
%

function [aap]=aas_addtask(aap,taskname,varargin)
v=varargin;
if (~isempty(v) && isstruct(v{end}))
    extraparameters=v{end};
    v{end}=[];
else
    extraparameters=[];
end;

if (~isempty(v))
    epiprefix=v{1};
    v{1}=[];
else
    epiprefix=[];
end;

if (~isempty(v))
    tobecompletedfirst=v{1};
    v{1}=[];
else
    tobecompletedfirst='[previous]';
end;

[aap index]=aas_addtaskparameters(aap,taskname);
mystage.index=index;
mystage.name=taskname;

if isfield(aap.tasksettings.(taskname)(index),'specialrequirements') && ...
        isfield(aap.tasksettings.(taskname)(index).specialrequirements,'initialisationmodule')
    % add as initialisation module
    aap.tasklist.initialisation.module=[aap.tasklist.initialisation.module mystage];
    
else % add as main module
    
    mystage.epiprefix=epiprefix;
    mystage.tobecompletedfirst=tobecompletedfirst;
    mystage.extraparameters=extraparameters;
    mystage.remotestream=[];
    
    if isempty(aap.tasklist.main.module) % start from empty list
        aap.tasklist.main(1).module(1)=mystage;
    elseif isempty(aap.tasklist.main.module(1).name) % replace empty task
        aap.tasklist.main(1).module(1)=mystage;
    else % append to existing tasks
        if ~isfield(mystage,'aliasfor'), mystage.aliasfor=[]; end % djm
        aap.tasklist.main.module=[aap.tasklist.main.module mystage];
        aap.aap_beforeuserchanges.tasklist.main.module=[aap.aap_beforeuserchanges.tasklist.main.module mystage];
    end
    
end
