% aa queue processor that uses Matlab's Parallel Computing Toolbox to run
% tasks.
%
% It differs from the Condor or qsub queue managers in that it uses a 
% pre-existing pool of matlab sessions, started with a prior command like 
% "matlabpool 8". As there's no need to launch a matlab session for each 
% processing module, latencies are reduced, which is important if you have
% lots of small jobs (i.e., like realtime)
%
% It also has more efficient dependency calculation, which rather than
% using done flags keeps track (within the current processing session).
% Done flags are stil read in at the beginning and written as a module
% finishes.

classdef aaq_matlab_pct<aaq
    properties
        maxretries=5;
        
        filestomonitor=[];
        filestomonitor_jobnum=[];
        compiledfile=[]; 
        matlab_pct_path=[];
        retrynum=[];
        jobnotrun=[];
        fatalerrors=false;
        jobcount=0;
        jobstatus=[];
        dep_names={};
        dep_done=[];
        depof={};
        depon={};
        workerstatus={};
    end
    methods
        function [obj]=aaq_matlab_pct(aap)
            global aaworker
            obj.aap=aap;
            
            obj.matlab_pct_path=fullfile(aaworker.parmpath,'matlab_pct');
            if (exist(obj.matlab_pct_path,'dir')==0)
                mkdir(obj.matlab_pct_path);
            end;
            
        end
        %%==============================
        % Add a task to the task queue
        % This adds to the parent class's method, by scanning the dependencies
        % of the new module.
        %
        function obj=addtask(obj,taskmask)
            obj=addtask@aaq(obj,taskmask);
            
            njob=length(obj.jobqueue);
            
            obj.retrynum(njob)=0;
            obj.jobnotrun(njob)=true;
            
            %% For better performance, make a list of all dependencies and refer to them by index
            % Removes need for endless done flag file checking
            obj.depof{njob}=[];
            obj.depon{njob}=[];
            for tbcfind=1:length(obj.jobqueue(njob).tobecompletedfirst)
                % Find or put this dependency in a list
                df=obj.jobqueue(njob).tobecompletedfirst{tbcfind};
                depind=find(strcmp(df,obj.dep_names));
                if isempty(depind)
                    obj.dep_names{end+1}=df;
                    depind=length(obj.dep_names);
                    obj.dep_done(depind)=exist(df,'file');
                end;
                
                % Only record dependencies that haven't already been
                % satisfied
                if ~obj.dep_done(depind)
                    obj.depon{njob}(end+1)=depind;
                    obj.depof{depind}(end+1)=njob;
                end;
            end;
            
        end
        
        
        
        % Run all jobs, using matlabs "single program, multiple data"
        % construct. This isn't really what we use it for - we actually set
        % pass messages between the modules to make a queue happen
        function [obj]=runall(obj,dontcloseexistingworkers,waitforalljobs)
            if ~isempty(obj.jobqueue)
                spmd
                    if labindex==1
                        obj.runall_spmd_master(waitforalljobs);
                    else
                        obj.runall_spmd_worker();
                    end;
                end;
            end;
        end;
        
        % Used by the master thread to pump any messages from workers
        function [obj,messagesreceived]=receiveWorkerMessages(obj)
            messagesreceived=false;
            
            while(labProbe)
                [data source tag]=labReceive;
                fprintf('Worker %d message %s which was %s\n',source,data{1},data{2});
                switch(data{1})
                    case 'status'
                        obj.workerstatus{source}=data{2};
                        fprintf('Worker %d status update %s\n',source,data{2});
                    case 'done_aa_doprocessing'
                        jobind=data{2};
                        obj.jobnotrun(jobind)=false;
                        aas_log(obj.aap,false,sprintf('PARALLEL (matlab_pct): Completed %s',obj.jobqueue(jobind).doneflag));
                        % Remove this as a dependency from all of the
                        % stages dependent on the stage that has just
                        % completed
                        for depoflist=1:length(obj.depof{jobind})
                            obj.depon{obj.depof{jobind}(depoflist)}(obj.depon{obj.depof{jobind}(depoflist)}==jobind)=[];
                        %    fprintf('Cleared dependency %d from %d, it has %d remaining\n',jobind,obj.depof{jobind}(depoflist),length(obj.depon{obj.depof{jobind}(depoflist)}));                          
                        end;
                    otherwise
                        fprintf('Unrecognized message from worker %s\n',data{1});

                end;
                messagesreceived=true;
            end;
        end;
        
        % Master - sends jobs that are ready to master
        function [obj]=runall_spmd_master(obj,waitforalljobs)
            obj.receiveWorkerMessages();
            
            itnum=0;
            
            fprintf('Now trying job execution\n');
            obj
            % Main job execution list
            while(any(obj.jobnotrun) || not(isempty(obj.filestomonitor)))
                itnum=itnum+1;
                for jobind=1:length(obj.jobqueue)
                    if obj.jobnotrun(jobind) && isempty(obj.depon{jobind})
                        fprintf('Job %d ready to go.\n',jobind);
                        % Job has no dependencies - good to go...
                        obj.jobcount=obj.jobcount+1;
                        job=obj.jobqueue(jobind);
                        obj.aap.acq_details.root=aas_getstudypath(obj.aap,job.k);
                        job.aap=obj.aap;
                        
                        fprintf('Waiting for a free worker\n');
                        % Find out who can run it, wait if no-one free
                        while(1)
                            workerwaiting=find(strcmp('waiting',obj.workerstatus));
                            if ~isempty(workerwaiting)
                                workerwaiting=workerwaiting(1);
                                obj.workerstatus{workerwaiting}='pending';
                                break;
                            end;
                            obj.receiveWorkerMessages();
                            pause(0.1);
                        end;
                        
                        
                        % Run the job
                        aas_log(obj.aap,false,sprintf('PARALLEL (matlab_pct): Submitting to worker %d job %s',workerwaiting,job.doneflag));
                        labSend({'aa_doprocessing',obj,job,jobind},workerwaiting);
                        
                        % Flag as submitted
                        obj.jobnotrun(jobind)=false;
                    end
                end
                
                if ~waitforalljobs
                    break;
                end;
            
                obj.receiveWorkerMessages();

                % Lets not overload the filesystem
                pause(0.1);
            end
            
            
            % Close all workers as they become free
            fprintf('Closing\n');
            while any(~strcmp('closed',obj.workerstatus))
                workerwaiting=find(strcmp('waiting',obj.workerstatus));
                for wwind=1:length(workerwaiting)
                    fprintf('Sending close to %d\n',workerwaiting(wwind));
                    labSend({'close'},workerwaiting(wwind));
                end;
                obj.receiveWorkerMessages();
                pause(0.1);
            end;
            
            
            fprintf('All done\n');
            if waitforalljobs == 1
                obj.emptyqueue;
            end
            
            if (obj.fatalerrors)
                aas_log(obj.aap,false,'PARALLEL (matlab_pct): Fatal errors executing jobs. The errors were:','Errors');
                for errind=1:length(errline)
                    aas_log(obj.aap,false,[' ' errline{errind}],'Errors');
                end;
                aas_log(obj.aap,true,'PARALLEL (matlab_pct): Terminating because of errors','Errors');
            end
        end
        
        
        
        %% Code for a worker - basically hangs around waiting to be told what 
        % to do, and then executes code when asked
        function [obj]=runall_spmd_worker(obj)
            alldone=false;
            while ~alldone
                % Say we're waiting for something to do
                labSend({'status','waiting'},1);
                
                % Wait to be told what that is
                while(~labProbe)
                    pause(0.1);
                end;
                
                
                [data source tag]=labReceive;
                switch(data{1})
                    case 'aa_doprocessing'
                        labSend({'status','aa_doprocessing'},1);
                        obj=data{2};
                        job=data{3};
                        jobind=data{4};
                        aa_doprocessing_onetask(obj.aap,job.task,job.k,job.indices);
                        labSend({'done_aa_doprocessing',jobind},1);
                    case 'close'
                        labSend({'status','closed'},1);
                        alldone=true;
                end;
            end;
        end;
        
    end;
    
end