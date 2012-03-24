% This is a message pump for the worker

function aaworker_poll(obj,event)
global aaworker;
numretries=1; % Now retries are handling by master - better at reallocating in case of memory errors etc.
stop(aaworker.polltimer);
mastertaskpth=fullfile(aaworker.parmpath,'pendingtask.mat');
taskpth=aas_propagatefrom(aaworker.master.hostname,mastertaskpth,'worker');
if (exist(taskpth))
    for i=1:30
        try
            task=load(taskpth);
            task=task.task;
            aaworker.aap=task.aap;
            break;
        catch
            % assume it was a file propagation delay or locking problem and so try again
            taskpth=aas_propagatefrom(aaworker.master.hostname,mastertaskpth,'worker');
            fprintf('waiting for propagation %d\n',i);
            pause(5.0);
        end;
    end;
    fprintf('AAWORKER: Found task %s here %s\n',task.name,taskpth);
    tic;
    switch (task.name)
        case 'doprocessing'
            for i=1:length(aaworker.aap.internal.taskqueue)

                % If at first you don't succeed, try, and try again
                for jj=1:numretries
                    try
                        mytask=aaworker.aap.internal.taskqueue(i);
                        % allow full path of module to be provided
                        [stagepath stagename]=fileparts(aaworker.aap.tasklist.main.module(mytask.k).name);
                        stagetag=aas_getstagetag(aaworker.aap,mytask.k);
                        try
                            aaworker.aap.tasklist.currenttask.settings=aaworker.aap.tasksettings.(stagename)(aaworker.aap.tasklist.main.module(mytask.k).index);
                        catch
                            aaworker.aap.tasklist.currenttask.settings=[];
                        end;

                        switch(mytask.domain)
                            % now run current stage
                            case 'study'
                                aas_log(aaworker.aap,0,sprintf('\nAAWORKER %s RUNNING: %s',stagetag,mytask.description));
                                [aaworker.aap,resp]=feval(aaworker.aap.tasklist.main.module(mytask.k).name,aaworker.aap,mytask.task);
                                writedoneflag(aaworker,mytask.doneflag);
                                aas_log(aaworker.aap,0,sprintf('\nAAWORKER %s COMPLETED',stagetag));
                            case 'subject'
                                aas_log(aaworker.aap,0,sprintf('\nAAWORKER %s RUNNING: %s for %s',stagetag,mytask.description,aas_getsubjname(aaworker.aap,mytask.i)));
                                [aaworker.aap,resp]=feval(aaworker.aap.tasklist.main.module(mytask.k).name,aaworker.aap,mytask.task,mytask.i);
                                writedoneflag(aaworker,mytask.doneflag);
                                aas_log(aaworker.aap,0,sprintf('\nAAWORKER %s COMPLETED',aaworker.aap.tasklist.main.module(mytask.k).name));
                            case 'session'
                                aas_log(aaworker.aap,0,sprintf('\nAAWORKER %s RUNNING: %s for %s ',aaworker.aap.tasklist.main.module(mytask.k).name,mytask.description,aas_getsessname(aaworker.aap,mytask.i,mytask.j)));
                                [aaworker.aap,resp]=feval(aaworker.aap.tasklist.main.module(mytask.k).name,aaworker.aap,mytask.task,mytask.i,mytask.j);
                                writedoneflag(aaworker,mytask.doneflag);
                                aas_log(aaworker.aap,0,sprintf('\nAAWORKER %s COMPLETED',aaworker.aap.tasklist.main.module(mytask.k).name));
                            case 'internal'
                                if ~mytask.i % need to setup loop variable
                                    [aaworker.aap,resp]=feval(aaworker.aap.tasklist.main.module(mytask.k).name,aaworker.aap,'parallelise');
                                    mkdir(mytask.doneflag);
                                    aas_propagateto(aaworker.master.hostname,mytask.doneflag);
                                    aas_log(aaworker.aap,0,sprintf('\nAAWORKER %s PARALLELISED',stagetag));
                                else
                                    aas_log(aaworker.aap,0,sprintf('\nAAWORKER %s RUNNING: %s',stagetag,mytask.description));
                                    [aaworker.aap,resp]=feval(aaworker.aap.tasklist.main.module(mytask.k).name,aaworker.aap,mytask.task,mytask.i);
                                    writedoneflag(aaworker,mytask.doneflag); % includes propagation
                                    aas_log(aaworker.aap,0,sprintf('\nAAWORKER %s COMPLETED',aaworker.aap.tasklist.main.module(mytask.k).name));
                                    doneflagpath=fileparts(mytask.doneflag);
                                    alldone=aas_checkinternaldomainprogress(doneflagpath);
                                    if alldone 
                                        writedoneflag(aaworker,fullfile(doneflagpath,'all.done'));
                                        aas_log(aaworker.aap,0,sprintf('\nALL %s SUBTASKS COMPLETED',aaworker.aap.tasklist.main.module(mytask.k).name));
                                    end
                                end
                        end;
                        jj=0;
                        % Flush diary to update
                        diary off
                        diary(aaworker.diaryname)
                        break;
                    catch;
                        % Flush diary to update
                        diary off
                        diary(aaworker.diaryname)
                        % make sure we're in the right directory
                        le=lasterror;
                        cd (aaworker.aap.internal.pwd);
                        aas_log(aaworker.aap,0,sprintf('AAWORKER WARNING: retrying as failed on attempt %d of %d with error:\n\t%s',jj,numretries,le.message));
                        for kk=1:length(le.stack)
                            aas_log(aaworker.aap,0,sprintf(' file %s name %s line %d',le.stack(kk).file,le.stack(kk).name,le.stack(kk).line));
                        end;
                    end;
                end; % next retry
                if (jj==numretries)
                    aas_log(aaworker.aap,1,sprintf('AAWORKER ERROR: %d retries did not suceed - see details above',numretries));
                end;
            end; % next task in queue
    end; % currently redundant
    try
        aaworker_setspmvisibility
        delete(taskpth);
        delete(mastertaskpth);
        timenow=now;
        boredfn=fullfile(aaworker.parmpath,'iambored.mat');
        save(boredfn,'timenow'); rehash
        aas_propagateto(aaworker.master.hostname,boredfn);
        aaworker.lastexcitement=clock;
    catch
        debugnow
        % Flush diary to update
        diary off
        diary(aaworker.diaryname)
    end;
else
    % Terminate worker if not used for 3 minutes
    if (etime(clock,aaworker.lastexcitement)>180)
        try
            aas_log(aaworker.aap,0,'There is nothing happening. I am out of here.\n');
            diary off
        catch
            % Flush diary to update
            diary off
            diary(aaworker.diaryname)
        end;
        % Now use recursive kill to finish self off, as this will get any
        % child processes
        aas_sackworker(aaworker.aap,aaworker.id);
        % bet you never get here
        exit;
    end;
end;
% Flush diary to update
diary off
diary(aaworker.diaryname)


start(aaworker.polltimer);

%%
function writedoneflag(aaworker,fn)
fid=fopen(fn,'w');
if (~fid) aas_log(aaworker.aap,1,['Error writing done flag ' fn]); end;
try fprintf(fid,'%f',toc);
catch
    keyboard
end
fclose(fid);
aas_propagateto(aaworker.master.hostname,fn);
return;



