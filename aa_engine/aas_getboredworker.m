function [aap workerid]=aas_getboredworker(aap,specialrequirements)
global aaparallel

if (~exist('specialrequirements','var'))
    specialrequirements={};
end;
if (isstruct(specialrequirements))
    specialrequirements=fieldnames(specialrequirements);
end;

% First check exisiting workers to see if any are bored
workerid=[];
for i=1:length(aaparallel.workerlist)
    % check if eligible for each job
    eligibleworker=true;
    for j=1:length(specialrequirements)
        if (~isfield(aaparallel.workerstatus.(sprintf('worker%d',aaparallel.workerlist(i))),'specialrequirements'))
            eligibleworker=false;
        else
            if ~isfield(aaparallel.workerstatus.(sprintf('worker%d',aaparallel.workerlist(i))).specialrequirements,(specialrequirements{j}))
                eligibleworker=false;
            end;
        end;
    end;
    % yes, got an eligible worker
    if (eligibleworker && strcmp(aaparallel.workerstatus.(sprintf('worker%d',aaparallel.workerlist(i))).status,'bored'))
        workerid=aaparallel.workerlist(i);
        aas_log(aap,false,sprintf('PARALLEL-MANAGE-WORKERS: got eligible worker (%d special requirements) worker%d',length(specialrequirements),workerid));
        aaparallel.workerstatus.(sprintf('worker%d',workerid)).status='joballocated';
        break;
    end;
end;

% If allowed, make a new worker

if (length(aaparallel.workerlist)<aaparallel.numberofworkers)
    workerid=aaparallel.nextworkernumber;
    aaparallel.nextworkernumber=aaparallel.nextworkernumber+1;
    aas_launchworker(aap,workerid,specialrequirements);
    % workers that don't offer special requirements should be used first,
    % and the special ones saved until they're needed. So add non-special
    % workers to the front of the worklist; special ones to the back
    if (isempty(specialrequirements))
        aaparallel.workerlist=[workerid aaparallel.workerlist];
    else
        aaparallel.workerlist=[aaparallel.workerlist workerid];
    end;
    pth=aaworker_getparmpath(aap,workerid);
    aaparallel.workerstatus.(sprintf('worker%d',workerid)).status='starting';
    aaparallel.workerstatus.(sprintf('worker%d',workerid)).statuschanged=clock;
    aaparallel.workerstatus.(sprintf('worker%d',workerid)).allocatedjobs=[];
    aaparallel.workerstatus.(sprintf('worker%d',workerid)).runningjob=[];
    for j=1:length(specialrequirements)
        aaparallel.workerstatus.(sprintf('worker%d',workerid)).specialrequirements.(specialrequirements{j})=1;
    end;
    aaparallel.workerstatus.(sprintf('worker%d',workerid)).status='joballocated';
end;
