%% Automatic analysis
% Update worker status
% Check the status of each of the workers and update their status fields.
% Also, if a worker with an allocated job becomes bored, it gets passed to
% them right away.
%

function aap=aas_updateworkerstatus(aap)
global aaparallel

aaparallel.workersummary.bored=0;
aaparallel.workersummary.joballocated=0;
aaparallel.workersummary.busy=0;
aaparallel.workersummary.starting=0;
aaparallel.workersummary.exited=0;
aaparallel.workersummary.unknown=0;
aaparallel.workersummary.error=0;
aaparallel.workersummary.total=length(aaparallel.workerlist);
i=1;

while(i<=length(aaparallel.workerlist))
    workerid=aaparallel.workerlist(i);
    numallocjobs=length(aaparallel.workerstatus.(sprintf('worker%d',workerid)).allocatedjobs);
    if (numallocjobs>1)
        aas_log(aap,false,sprintf('PARALLEL: unexpected, %d jobs allocated to a single worker %d',numallocjobs,workerid));
    end;
    oldstatus=aaparallel.workerstatus.(sprintf('worker%d',workerid)).status;
    pth=aaworker_getparmpath(aap,aaparallel.workerlist(i),true);
    if (~exist(pth,'file'))
        status='exited';
        aas_log(aap,false,sprintf('PARALLEL: worker %d exited at %s',workerid,datestr(now,14)));
        aaparallel.workerlist(i)=[];
        i=i-1;
    elseif (~exist(fullfile(pth,'aaworker.xml'),'file'))
        status='starting';
    elseif (exist(fullfile(pth,'errorflag.txt'),'file'))
        status='error';
        try
            cmd=['tail -n8 ' fullfile(pth,'log.txt')];
            [s w]=aas_shell(cmd);
            aaparallel.workerstatus.(sprintf('worker%d',workerid)).errorlog=w;
        catch
        end;
    elseif (exist(fullfile(pth,'iambored.mat'),'file'))
        if (~isfield(aaparallel.workerstatus.(sprintf('worker%d',workerid)),'workerdetails'))
            aaparallel.workerstatus.(sprintf('worker%d',workerid)).workerdetails=xml_read(fullfile(pth,'aaworker.xml'));
        end;
        if ~isempty(aaparallel.workerstatus.(sprintf('worker%d',workerid)).runningjob); %djm
            task=aaparallel.workerstatus.(sprintf('worker%d',workerid)).runningjob.aap.internal.taskqueue;
            if strcmp(task.domain,'internal') && ~task.i
                % parallelised and ready to go, so need to add new jobs
                aap=addnewjobs(aap,task);  
            end
        end
        aaparallel.workerstatus.(sprintf('worker%d',workerid)).runningjob=[];      
        if (length(aaparallel.workerstatus.(sprintf('worker%d',workerid)).allocatedjobs)>0)
            % Now communicate with the worker
            try
                delete(fullfile(pth,'iambored.mat'));
                task=aaparallel.workerstatus.(sprintf('worker%d',workerid)).allocatedjobs{1};
                save(fullfile(pth,'pendingtask'),'task');
                stagename=task.aap.internal.taskqueue.stagename;
                description=task.aap.internal.taskqueue.description;
                aas_log(aap,false,sprintf('PARALLEL: Worker %d sent %s: %s',workerid,stagename,description));

                aaparallel.workerstatus.(sprintf('worker%d',workerid)).runningjob= aaparallel.workerstatus.(sprintf('worker%d',workerid)).allocatedjobs{1};
                aaparallel.workerstatus.(sprintf('worker%d',workerid)).allocatedjobs(1)=[];
                status='busy';
            catch
                status='unknown';
            end;
        else
            status='bored';
        end;
    elseif (exist(fullfile(pth,'pendingtask.mat'),'file'))
        status='busy';
    else
        status='unknown';
    end;

    aaparallel.workerstatus.(sprintf('worker%d',workerid)).status=status;

    if (~strcmp(oldstatus,status))
        aaparallel.workerstatus.(sprintf('worker%d',workerid)).statuschanged=clock;
    end;

    aaparallel.workersummary.total=length(aaparallel.workerlist);

    aaparallel.workersummary.(status)=aaparallel.workersummary.(status)+1;

    i=i+1;
end;

%% add new jobs for newly parallelised domain [djm]
function aap=addnewjobs(aap,task)

doneflagdir=aas_doneflag_getpath(aap,task.k);

alldone=aas_checkinternaldomainprogress(doneflagdir);

if alldone
    alldoneflag=fullfile(doneflagdir,'all.done');
    fid=fopen(alldoneflag,'w');
    if (~fid); aas_log(aap,1,['Error writing done flag ' alldoneflag]); end;
    fprintf(fid,'%s',datestr(now));
    fclose(fid);
    aas_log(aap,0,sprintf('- all subjobs completed previously for: %s',task.stagename));
else
    % load joblist
    joblist=fullfile(aap.acq_details.root, ...
        sprintf('%s_%g.parallel.mat',task.stagename,aap.tasklist.main.module(task.k).index));
    tries=0;
    while tries<15
        try 
          rehash; 
          load(joblist);
          break
        catch
            % file system too slow? Have seen it take >10 seconds for the file to become visible
            pause(2); tries=tries+1;
        end
    end

    for x=1:length(loopvar)
        if ischar(loopvar{x})
            doneflag=fullfile(doneflagdir,[loopvar{x} '.done']); % file in directory
            desc=loopvar{x};
        elseif isstruct(loopvar{x})
            try c={loopvar{x}.id};
            catch c=strcat(struct2cell(loopvar{x}),'.'); % Convert structure to cell
            end
            doneflag=fullfile(doneflagdir,[c{:} 'done']); % file in directory
            desc=[c{:}];
        end
        if (length(dir(doneflag)))
            aas_log(aap,0,sprintf('- completed previously: %s for %s',task.stagename,desc));
        else
            % now queue current stage
            taskmask=task;
            taskmask.i=x;
            taskmask.tobecompletedfirst=[];
            taskmask.doneflag=doneflag;
            taskmask.description=sprintf('%s for %s',task.description,desc);
            aap.internal.taskqueue=[aap.internal.taskqueue,taskmask];
            aas_log(aap,0,sprintf('MODULE %s QUEUED: for %s',task.stagename,desc));
        end;
    end % next subjob
    
end % create all done flag, or queue subjobs?

return
