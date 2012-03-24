classdef aaq_localparallel<aaq
    methods
        function [obj]=aaq_localparallel(aap)
            obj.aap=aap;
        end;

%% ==============================
        % Run all tasks on the queue
        function [obj]=runall(obj,dontcloseexistingworkers)
            global aaparallel
            
            njobs=length(obj.jobqueue);
            
            % FIRST SECTION CHECKS IF JOBS ARE ALREADY RUNNINIG
            if (strcmp(dontcloseexistingworkers,'cont') || strcmp(dontcloseexistingworkers,'continue'))
                % update status
                aas_updateworkerstatus(obj.aap);
                
                % First of all, find out which jobs are ready to run
                for i=1:njobs
                    % Find out whether this job is ready to be allocated
                    for j=1:length(obj.jobqueue(i).tobecompletedfirst)
                        if (~exist(obj.jobqueue(i).tobecompletedfirst{j},'file'))
                            obj.jobqueue(i).readytorun=false;
                        end;
                    end;
                end;
                
                % For all jobs ready to run, make a task description
                taskdesc=cell(njobs,1);
                for i=1:njobs
                    if (obj.jobs(i).readytorun)
                        taskdesc{i}=obj.jobqueue(i).getjobdescription(obj.aap,i);
                        taskdesc{i}.taskqueueposition='notimportant';
                    end;
                end;
                
                % Go through each of the existing workers, and find out if they're already doing
                % one of the tasks we're about to try to allocate
                i=1;
                while(i<=length(aaparallel.workerlist))
                    workerid=aaparallel.workerlist(i);
                    %                    oldstatus=aaparallel.workerstatus.(sprintf('worker%d',workerid)).status;
                    pth=aaworker_getparmpath(obj.aap,workerid,true);
                    workeruseful=false;
                    try
                        task=load(fullfile(pth,'pendingtask'),'task');
                        task=task.task;
                        task.taskqueueposition='notimportant';
                        for j=1:njobs
                            if (obj.jobs(j).readytorun)
                                if (isequal(taskdesc{j},task))
                                    mytaskdesc=task.obj.jobqueue.description;
                                    aas_log(obj.aap,false,sprintf('PARALLEL: Will not sack worker %d as already on useful task %s',workerid,mytaskdesc));
                                    obj.jobs(j).joballocated=true;
                                    obj.jobs(j).lastrunat=now;
                                    obj.jobs(j).timesattempted=obj.jobs(j).timesattempted+1;
                                    workeruseful=true;
                                    break;
                                end;
                            end;
                        end;
                    catch
                    end;
                    if (~workeruseful)
                        aas_sackworker(obj.aap,workerid);
                    else
                        i=i+1;
                    end;
                end;
            end
            
            
            %% NOW ACTUALLY RUN THEM, LAUNCHING NEW WORKERS IF POSSIBLE
            pendingjobs=99999;
            awaitingprevious=99999;
            olddesc=0;
            seriouserrors=[];
            retrydelays=aaparallel.retrydelays/24/3600; % convert to days
            numtries=length(retrydelays);
            
            while (~isempty([awaitingprevious pendingjobs])) || ~isempty(aaparallel.workerlist)
                awaitingprevious=[];
                pendingjobs=[];
                numallocated=[];
                awaitingretry=[];
                manage_workers_log=false;
                
                for i=1:njobs
                    if (~obj.jobs(j).joballocated && ~ismember(i,seriouserrors))
                        % Find out whether this job is ready to be allocated by
                        % checking dependencies (done_ flags)
                        readytorun=true;
                        for j=1:length(obj.jobqueue(i).tobecompletedfirst)
                            if (~exist(obj.jobqueue(i).tobecompletedfirst{j},'file'))
                                readytorun=false;
                            end;
                        end;
                        
                        % If retrying, has an elegant delay elapsed?
                        notyetretrytime=false;
                        if ((timesattempted(i)>0) && (timesattempted(i)<=numtries))
                            if (now<(lastrunat(i)+retrydelays(timesattempted(i))))
                                notyetretrytime=true;
                                readytorun=false;
                            end;
                        end;
                        
                        % If so, try to allocate it
                        if (readytorun)
                            pendingjobs=[pendingjobs i];
                            if (timesattempted(i)==(1+numtries))
                                aas_log(obj.aap,false,sprintf('***SERIOUS ERROR, Job failed %d times, giving up on it.',numtries));
                                seriouserrors=[seriouserrors i];
                            elseif (timesattempted(i)<=numtries)
                                if (timesattempted(i)>1) % djm: force to use high memory machine
                                    [obj obj.jobs(j).joballocated]=obj.allocate(i,'tryhighmemorymachine');
                                else % respect highmemory setting
                                    [obj obj.jobs(j).joballocated]=obj.allocate(i);
                                end
                                if (obj.jobs(j).joballocated)
                                    obj.jobs(i).timesattempted=obj.jobs(i).timesattempted+1;
                                    if (obj.jobs(i).timesattempted>1)
                                        aas_log(obj.aap,false,sprintf('PARALLEL: attempt %d, previous attempt was %d s ago.',timesattempted(i),round(3600*24*(now-lastrunat(i)))));
                                    end;
                                    obj.jobs(i).lastrunat=now;
                                else
                                    if (~manage_workers_log)
                                        aas_log(obj.aap,false,[sprintf('PARALLEL-MANAGE-WORKERS: workforce at capacity (%d/%d), cannot ask for more.',length(aaparallel.workerlist),aaparallel.numberofworkers)]);
                                        manage_workers_log=true;
                                    end;
                                end;
                            end;
                        else
                            if (notyetretrytime)
                                awaitingretry=[awaitingretry i];
                            else
                                awaitingprevious=[awaitingprevious i];
                            end;
                        end;
                    else
                        numallocated=[numallocated i];
                    end;
                end;
                
                % Now, actually provide the jobs that have been allocated to the
                % workers
                obj.aap=aas_updateworkerstatus(obj.aap);
                
                i=1;
                while (i<=length(aaparallel.workerlist))
                    sackthisworker=false;
                    workerid=aaparallel.workerlist(i);
                    % Now lets see decide upon any timeouts
                    sincestatuschanged=etime(clock,aaparallel.workerstatus.(sprintf('worker%d',workerid)).statuschanged);
                    try
                        timeout=60*obj.aap.timeouts.(aaparallel.workerstatus.(sprintf('worker%d',workerid)).status);
                        % Its been too long: give up and mark job as not allocated
                        if (timeout<sincestatuschanged)
                            sackthisworker=true;
                            sackedwhy=sprintf('Giving up on %s phase',aaparallel.workerstatus.(sprintf('worker%d',workerid)).status);
                            
                        end;
                    catch
                    end;
                    
                    % Is there an error?
                    if (strcmp(aaparallel.workerstatus.(sprintf('worker%d',workerid)).status,'error'))
                        sackthisworker=true;
                        sackedwhy='Error signal received';
                    end;
                    
                    % Is the process still alive?
                    try
                        pid=aaparallel.workerstatus.(sprintf('worker%d',workerid)).workerdetails.pid(1);
                        cmd=['ssh -x ' aaparallel.workerstatus.(sprintf('worker%d',workerid)).workerdetails.hostname sprintf(' ps -e | grep %d',pid)];
                        [s w]=aas_shell(cmd);
                        w=strtrim(w);
                        if (isempty(w))
                            sackthisworker=true;
                            sackedwhy='Process no longer running';
                        end;
                        if (~isempty(findstr(w,'Interrupt')))
                            aas_log(obj.aap,true, 'Ctrl-C pressed, exiting');
                        end;
                    catch
                    end;
                    
                    if (strcmp(aaparallel.workerstatus.(sprintf('worker%d',workerid)).status,'error'))
                        sackthisworker=true;
                        try
                            myerror=aaparallel.workerstatus.(sprintf('worker%d',workerid)).errorlog;
                        catch;
                            myerror='';
                        end;
                        sackedwhy=sprintf('Error signal received, which was\n>>>>>>\n%s\n<<<< ',myerror);
                    end;
                    
                    % Has it exited?
                    if (strcmp(aaparallel.workerstatus.(sprintf('worker%d',workerid)).status,'exited'))
                        sackthisworker=true;
                        sackedwhy='Worker exited unexpectedly';
                    end;
                    
                    if (sackthisworker)
                        aas_log(obj.aap,false,sprintf('PARALLEL: %s - sacking worker %d and reallocating any jobs',sackedwhy,workerid));
                        
                        aaparallel.workerstatus.(sprintf('worker%d',workerid)).status='sacked';
                        obj.aap=aas_sackworker(obj.aap,workerid);
                        alljobs=[aaparallel.workerstatus.(sprintf('worker%d',workerid)).allocatedjobs aaparallel.workerstatus.(sprintf('worker%d',workerid)).runningjob];
                        if (~~isempty(alljobs))
                            aas_log(obj.aap,false,'PARALLEL: no jobs to reallocate');
                        else
                            reamsg='PARALLEL: reallocating';
                            for j=1:length(alljobs)
                                jobnumber=alljobs{j}.taskqueueposition;
                                newerror=[];
                                newerror.sackedwhy=sackedwhy;
                                newerror.when=clock;
                                errorlist{jobnumber}=[errorlist{jobnumber} newerror];
                                joberrordir=fullfile(obj.aap.acq_details.root,'joberrors');
                                aas_makedir(obj.aap,joberrordir);
                                joberrorfn=fullfile(joberrordir,sprintf('J%d_%d.html',jobnumber,length(errorlist{jobnumber})));
                                ht=[sprintf('<html><body><h1>Error for J%d attempt %d</h1>', jobnumber,length(errorlist{jobnumber}))];
                                ht=[ht '<h2>' obj.jobqueue(jobnumber).description '</h2>'];
                                sackedwhy=strrep(sackedwhy,'\n','<br>');
                                ht=[ht sackedwhy];
                                ht=[ht '</body></html>'];
                                fid=fopen(joberrorfn,'w');
                                fprintf(fid,'%s',ht);
                                fclose(fid);
                                reamsg=sprintf('%s\tJ%d',reamsg,alljobs{j}.taskqueueposition);
                                joballocated(alljobs{j}.taskqueueposition)=false;
                            end;
                            aas_log(obj.aap,false,reamsg);
                        end;
                    else
                        i=i+1;
                    end;
                end;
                workersummary=[]; workersummary.bored=[]; workersummary.starting=[]; workersummary.busy=[]; workersummary.joballocated=[];
                for i=1:length(aaparallel.workerlist)
                    workerid=aaparallel.workerlist(i);
                    status=aaparallel.workerstatus.(sprintf('worker%d',workerid)).status;
                    if (~isfield(workersummary,status))
                        workersummary.(status)=workerid;
                    else
                        workersummary.(status)=[workersummary.(status) workerid];
                    end;
                end;
                desc=[length(awaitingprevious) length(awaitingretry) length(pendingjobs)-length(seriouserrors) length(numallocated) length(seriouserrors) length(workersummary.joballocated) length(workersummary.starting) length(workersummary.busy) length(workersummary.bored) length(aaparallel.workerlist)];
                if (any(desc-olddesc))
                    ht='<HTML><HEAD><title>aa parallel job and worker status</title>';
                    ht=[ht '<meta http-equiv="refresh" content="10">'];
                    ht=[ht '</head>'];
                    ht=[ht '<body><h1>aa parallel status ' sprintf('%d-%d-%d %d:%d:%d',round(clock)) '</h1>'];
                    ht=[ht '<table><tr>'];
                    
                    ht=[ht '<td valign="top"><h2>Jobs</h2>'];
                    ht=[ht '<h3>Status</h3>'];
                    ht=[ht '<table><tr>'];
                    ht=[ht  '<td valign="top">' aas_getjobdescription(obj.aap,'awaiting previous stages',awaitingprevious,errorlist) '</td>'];
                    ht=[ht  '<td valign="top">' aas_getjobdescription(obj.aap,'awaiting retry time',setdiff(awaitingretry,seriouserrors),errorlist) '</td>'];
                    ht=[ht  '<td valign="top">' aas_getjobdescription(obj.aap,'unallocated',setdiff(pendingjobs,seriouserrors),errorlist) '</td>'];
                    ht=[ht  '<td valign="top">' aas_getjobdescription(obj.aap,'allocated or complete',numallocated,errorlist) '</td>'];
                    ht=[ht  '<td valign="top">' aas_getjobdescription(obj.aap,'serious errors?',seriouserrors,errorlist) '</td>'];
                    ht=[ht '</tr></table></td>'];
                    ht=[ht '<td valign="top"><h2>Workers</h2>'];
                    ht=[ht '<table><tr><td>'];
                    ht=[ht '<h3>Status</h3>'];
                    ht=[ht '<table><tr>'];
                    ht=[ht '<td valign="top">' aas_getworkerdescription(obj.aap,'allocated',workersummary.joballocated)  '</td>'];
                    ht=[ht '<td valign="top">' aas_getworkerdescription(obj.aap,'starting',workersummary.starting)  '</td>'];
                    ht=[ht '<td valign="top">' aas_getworkerdescription(obj.aap,'busy',workersummary.busy) '</td>' ];
                    ht=[ht '<td valign="top">' aas_getworkerdescription(obj.aap,'bored',workersummary.bored) '</td>'];
                    ht=[ht '</tr></table>'];
                    
                    ht=[ht '</td></tr><tr><td>'];
                    ht=[ht '<h3>Recent management</h3>'];
                    if (isfield(aaparallel,'manage_workers_log'))
                        ht=[ht '<table>'];
                        for i=1:min(10,length(aaparallel.manage_workers_log))
                            ht=[ht '<tr>'];
                            ht=[ht '<td>' sprintf('%d-%d-%d %d:%d:%d',round(aaparallel.manage_workers_log((end+1-i)).when)) '</td>'];
                            ht=[ht '<td>' aaparallel.manage_workers_log((end+1-i)).msg '</td>'];
                            ht=[ht '</tr>'];
                        end;
                        ht=[ht '</table>'];
                    end;
                    ht=[ht '</td></tr></table>'];
                    ht=[ht '</td></tr></table>'];
                    ht=[ht '</body></html>'];
                    statusfn=fullfile(obj.aap.acq_details.root,'aa_parallel_status.html');
                    fid=fopen(statusfn,'w');
                    
                    %%% djm: keep trying for a minute
                    tries=0;
                    while fid==-1 && tries<30
                        pause(2)
                        fid=fopen(statusfn,'w');
                        tries=tries+1;
                    end
                    
                    fprintf(fid,'%s',ht);
                    fclose(fid);
                    olddesc=desc;
                end;
                
                pause(0.1);
            end;
            aas_log(obj.aap,false,'PARALLEL all jobs submitted to workers, waiting for completion');
            
            obj.emptyqueue;
        end;
        
        
        end;
        end
        