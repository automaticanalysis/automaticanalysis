% New version Mar 2008
% Doesn't wait for worker to become ready before returning - speeds up
% launching, but requires more careful checking by calling routine

function [aap workerid]=aas_getboredworker(aap)
global aaparallel

% First check exisiting workers to see if any are bored
workerid=[];

i=1;
while(i<=length(aaparallel.workerlist))
    % Check not already allocated to some other job that is waiting for it
    alreadyallocated=false;
    for i=1:length(aap.internal.taskqueue)
        if (aap.internal.taskqueue(i).jobstatus==2 && aap.internal.taskqueue(i).workerid==aaparallel.workerlist(i))
            alreadyallocated=true;
        end;
    end;
    if (~alreadyallocated)
        pth=aaworker_getparmpath(aap,aaparallel.workerlist(i),true);
        if (exist(pth,'dir'))
            if (~exist(fullfile(pth,'pendingtask.mat'),'file') && exist(fullfile(pth,'iambored.mat'),'file'))
                workerid=aaparallel.workerlist(i);
                break;
            end;
        else
            aas_log(aap,0,sprintf('PARALLEL worker %d exited.',aaparallel.workerlist(i)));
            aaparallel.workerlist(i)=[];
        end;
    end;
    i=i+1;
end;

if ((length(workerid)==0) & (length(aaparallel.workerlist)<aaparallel.numberofworkers))
    workerid=aaparallel.nextworkernumber;
    aaparallel.nextworkernumber=aaparallel.nextworkernumber+1;
    aas_launchworker(aap,workerid);
    aaparallel.workerlist=[aaparallel.workerlist workerid];
end;
