classdef aaq_matlab_pct < aaq
    % aaq_matlab_pct < aaq
    % aa queue processor that uses Matlab's spmd construct to run tasks.
    %
    % It differs from the Condor or qsub queue managers in that it uses a
    % pre-existing pool of matlab sessions, started with a prior command
    % like "parpool open 8". As there's no need to launch a matlab session
    % for each processing module, latencies are reduced, which is important
    % if you have lots of small jobs (i.e., like realtime)
    %
    % It also has more efficient dependency calculation, which rather than
    % using done flags keeps track (within the current processing session).
    % Done flags are still read in at the beginning and written as a module
    % finishes.
    %
    %
    % The sequence of method/function calls triggering the execution of a
    % module's code is
    %
    %   runall > runall_spmd_master | runall_spmd_worker > aa_doprocessing_onetask
    % 
    %
    % aaq_matlab_pct Properties (*inherited):
    %   doHandlePool - logical, true indicating that this class handles
    %                  the parpool
    %   benchmark    - struct, container for performance statistics
    %   matlab_pct_path - char arr, subdirectory for storing job diaries
    %   jobcount     - double, counter for jobs in pipeline
    %   *aap         - struct, the 'main' aa struct govering analysis
    %   *isOpen      - logical, flag indicating whether taskqueue is open
    %   *fatalerrors - logical, flag indicating fatal error
    %   *jobqueue    - struct, composed of 'taskmasks'
    % 
    % aaq_matlab_pct Methods (*inherited):
    %   addtask             - Call superclass' addtask and additionally 
    %                         scan for dependencies
    %   runall              - Run all jobs/tasks on the queue using spmd
    %   runall_spmd_master  - Send data to workers for processing (will be
    %                         run by worker with labindex 1, the 'master')
    %   runall_spmd_worker  - Wait for signal from master and either call 
    %                         aa_doprocessing_onetask or signal a closed
    %                         state (will be run by non-master workers)
    %   receiveWorkerMessages - Receive worker message and process it
    %   processWorkerMessage - Set workerstatus property or, if message is
    %                          'done_aa_doprocessing', perform wrap-up
    %   waitforrealtime     - <seemingly not used>
    %   close               - Close parpool and call superclass' close
    %   *save               - save self to file
    %   *emptyqueue         - clear the task queue (set jobqueue to [])
    %   *allocate           - allocate a job from the task queue to a worker
    %   *getjobdescription  - return struct task

    
    % ==================== 'Hidden help' for developers ===================
    % The worker with labindex 1 is defined as the master; communication
    % between master and the other workers is accomplished via labSend,
    % labReceive, and labProbe. 
    %
    % The master sends instructions and data to the other workers when they
    % are available:
    %   labSend({'aa_doprocessing',obj,job,...
    %   labSend({'close'},...
    %
    % The other workers always report back to the master (note the
    % labindex), in most cases reporting their status:
    %   labSend({'status','waiting'},1);
    %   labSend({'status','aa_doprocessing'},1);
    %   labSend({'status','closed'},1);
    % However, once a job is successfully finished, instead of the string
    % 'status' the worker in question sends 'done_aa_doprocessing', the job
    % index and a stream to return (which may be empty):
    %   labSend({'done_aa_doprocessing',jobind,streamcache_toreturn},1);
    % 
    % Hidden property 'workerstatus' is central for this communication in
    % that it indicates the current status of each worker ('pending',
    % 'waiting', 'waitforrealtime', 'idle', 'closed').
    % 
    % Above-mentioned more efficient dependency calculation is done in
    % method addtask, using hidden properties isJobReadyToGo, isJobNotRun,
    % jobDoneFlag, isJobDoneFlag, jobDepOn, jobDepOf, numJobDepOn.
    
    properties
        doHandlePool    = false; % logical, true indicating that class handles parpool
        benchmark       = struct('receiveWorkerMessages',0,...
                            'submitJobs',0,...
                            'waitForWorkers',0,...
                            'labsend',0,...
                            'jobreadystart',0); % struct, container for performance statistics
        matlab_pct_path = []; % subdirectory for storing job diaries (<aaworker.parmpath>/matlab_pct)
        jobcount        = 0;  % counter for jobs in pipeline
    end
    
    properties (Hidden)
        jobStudyPaths           = {} % cell array, actualised studypath (for each job)
        isJobReadyToGo          = [] % logical array, true indicating ready to be processed (for each job)
        isJobNotRun             = [] % logical array, true indicating job has not run (for each job)
        jobDoneFlag             = {} % cell array, copy of taskmask.doneflag (for each job)
        isJobDoneFlag           = [] % logical array, true if doneflag file exists (for each job)
        jobDepOn                = {} % cell array, jobs the actual job depends on (for each job)
        jobDepOf                = {} % cell array, jobs depending on the actual job (for each job)
        numJobDepOn             = [] % double array, no. of jobs the actual job depends on (for each job)

        aaworker                = [] % copy of global aaworker struct
        aaparallel              = [] % copy of global aaparallel struct
        workerstatus            = {} % cell array of chars indicating worker status (for each worker)
        
        % real-time
        realtime_deps           = [] % struct, dependencies of real-time jobs
        % Torque-only
        initialSubmitArguments  = '' % additional arguments to use when submitting jobs
    end
    
    methods
        function [obj]=aaq_matlab_pct(aap)
            % Constructor: setting up engine
            global aaparallel
            global aaworker
            obj.aap=aap;
            obj.aaworker = aaworker;
            obj.aaparallel = aaparallel;
            
            obj.matlab_pct_path=fullfile(obj.aaworker.parmpath,'matlab_pct');
            if ~exist(obj.matlab_pct_path,'dir')
                aas_makedir(obj.aap,obj.matlab_pct_path);
            end
            
            if ~aas_matlabpool('isopen')
                if ~isempty(aap.directory_conventions.poolprofile)
                    profiles = parallel.clusterProfiles;
                    if ~any(strcmp(profiles,aap.directory_conventions.poolprofile))
                        ppfname = which(spm_file(aap.directory_conventions.poolprofile,'ext','.settings'));
                        if isempty(ppfname)
                            aas_log(aap,true,sprintf('ERROR: settings for pool profile %s not found!',aap.directory_conventions.poolprofile));
                        else                            
                            P=parcluster(parallel.importProfile(ppfname));
                        end
                    else
                        aas_log(aap,false,sprintf('INFO: pool profile %s found',aap.directory_conventions.poolprofile));
                        P=parcluster(aap.directory_conventions.poolprofile);
                    end
                    switch class(P)
                        case 'parallel.cluster.Torque'
                            aas_log(aap,false,'INFO: Torque engine is detected');
                            P.ResourceTemplate = sprintf('-l nodes=^N^,mem=%dGB,walltime=%d:00:00', obj.aaparallel.memory,obj.aaparallel.walltime);
                            if any(strcmp({aap.tasklist.main.module.name},'aamod_meg_maxfilt')) && ... % maxfilt module detected
                                    ~isempty(aap.directory_conventions.neuromagdir) % neuromag specified
                                obj.initialSubmitArguments = ' -W x=\"NODESET:ONEOF:FEATURES:MAXFILTER\"';
                            end
                            P.SubmitArguments = strcat(P.SubmitArguments,obj.initialSubmitArguments);
                        case 'parallel.cluster.Generic'
                            aas_log(aap,false,'INFO: Generic engine is detected');
                            P.CommunicatingSubmitFcn = obj.SetArg(P.CommunicatingSubmitFcn,'walltime',obj.aaparallel.walltime);
                            P.CommunicatingSubmitFcn = obj.SetArg(P.CommunicatingSubmitFcn,'memory',obj.aaparallel.memory);                            
                        case 'Local'
                            aas_log(obj.aap,false,'INFO: Local engine is detected');
                    end
                else
                    P = parcluster('local');
                end
                P.NumWorkers = obj.aaparallel.numberofworkers;
                P.JobStorageLocation = obj.aaworker.parmpath;
                C = aas_matlabpool(P,P.NumWorkers);
                if ~isempty(C)
                    C.IdleTimeout = obj.aaparallel.walltime*60; 
                end
                obj.doHandlePool = true;                
            end
        end

        
        function close(obj)
            % Closes matlabpool and invokes superclass' close function.
            if obj.doHandlePool
                aas_matlabpool('close'); 
            end
            close@aaq(obj);
        end
        
        
        function obj=addtask(obj,taskmask)
            % Add a task to the task queue.
            % Adds to the parent class's method by scanning the
            % dependencies of the new module.
            obj=addtask@aaq(obj,taskmask);
            % module_index is needed for dealing with real-time input only
            module_index=obj.aap.tasklist.main.module(taskmask.k).index;
            job_ix=length(obj.jobqueue);
            obj.isJobNotRun(job_ix)=true;
            
            % For better performance, make a list of all dependencies and
            % refer to them by index 
            % Removes need for endless done flag file checking
            obj.jobDepOf{job_ix}=[];
            obj.jobDepOn{job_ix}=[];
            obj.numJobDepOn(job_ix)=0;
            
            % Store this task as a future dependency of others
            obj.jobDoneFlag{job_ix}=taskmask.doneflag;
            obj.isJobDoneFlag(job_ix)=exist(taskmask.doneflag,'file');
            obj.jobStudyPaths{job_ix}=aas_getstudypath(obj.aap,taskmask.k);
            
            % Does this stage need realtime input?
            if isfield(obj.aap.schema.tasksettings.(taskmask.stagename)(module_index).ATTRIBUTE,'waitforrealtime_singlefile') && ~isempty(obj.aap.schema.tasksettings.(taskmask.stagename)(module_index).ATTRIBUTE.waitforrealtime_singlefile)
                [dcmfield, dcmfilter]=strtok(obj.aap.schema.tasksettings.(taskmask.stagename)(module_index).ATTRIBUTE.waitforrealtime_singlefile,'=');
                newrtd=struct('dcmfield',dcmfield,'dcmfilter',strtrim(dcmfilter(2:end)),'njob',job_ix,'eventtype','singlefile','satisfied',false);
                if isempty(obj.realtime_deps)
                    obj.realtime_deps=newrtd;
                else
                    obj.realtime_deps(end+1)=newrtd;
                end
                obj.jobDepOn{job_ix}=-length(obj.realtime_deps); % minus indicates realtime dependency
                obj.numJobDepOn(job_ix)=1;
            end
            
            % Go through each dependency of this task and find the index of
            % that
            for tbcfind=1:length(obj.jobqueue(job_ix).tobecompletedfirst)
                % Find or put this dependency in a list
                df=obj.jobqueue(job_ix).tobecompletedfirst{tbcfind};
                depind=find(strcmp(df,obj.jobDoneFlag));
                if isempty(depind)
                    aas_log(aap,true,sprintf('Dependency %s of %s not found.',df,taskmask.doneflag));
                end
                
                % Only record dependencies that haven't already been
                % satisfied
                if ~obj.isJobDoneFlag(depind)
                    obj.jobDepOn{job_ix}(end+1)=depind;
                    obj.numJobDepOn(job_ix)=obj.numJobDepOn(job_ix)+1;
                    obj.jobDepOf{depind}(end+1)=job_ix;
                end
            end
            
            % Add it to the isJobReadyToGo queue if appropriate
            if obj.numJobDepOn(job_ix)==0
                obj.isJobReadyToGo=[obj.isJobReadyToGo job_ix];
            end
        end

        
        function [obj]=runall(obj, ~, waitforalljobs)
            % Run all jobs, using spmd.
            % This isn't really what we use it for - we actually set pass
            % messages between the modules to make a queue happen.
            % 'waitforalljobs' is a logical, true indicating that the job
            % queue is fully built, false otherwise.
            
            % This method only acts when the queue is fully built
            if waitforalljobs
                if ~isempty(obj.jobqueue)
                    % Copy current value of global variable aacache to
                    % obj.aaworker (necessary because workers invoked by
                    % spmd have 'local copies of global variables',
                    % initiated with [])
                    global aacache
                    obj.aaworker.aacache = aacache;
                    spmd
                        if labindex==1
                            obj.runall_spmd_master(waitforalljobs);
                        elseif labindex==2 && obj.aap.options.realtime.scanforincomingmeta
                            obj.waitforrealtime();
                        else
                            obj.runall_spmd_worker();
                        end
                    end
                end
            end
        end
        
        
        function [obj,messagesreceived]=receiveWorkerMessages(obj)
            % Pull data from workers and process them.
            % This method is used exclusively by the master worker.
            tic
            messagesreceived=false;
            while(labProbe)
                [workerData, workerIndex, ~]=labReceive;
                obj.processWorkerMessage(workerData,workerIndex);
                messagesreceived=true;
            end
            obj.benchmark.receiveWorkerMessages=obj.benchmark.receiveWorkerMessages+toc;
        end
        
        
        function [obj]=processWorkerMessage(obj,workerData,workerIndex)
            % Execute action implied by message from non-master worker.
            % Input arg 'workerData' is the data sent from any of the
            % non-master workers, see comments in the code.

            % the first element's value signals the action to be triggered:
            % if it is 'status', do nothing but copy the second element
            % to the appropriate position in obj.workerstatus; if it is
            % 'done_aa_doprocessing', 
            firstOutputArg = workerData{1};
            
            switch(firstOutputArg)
                case 'status'
                    % workerData{2} is a string indicating the worker's status,
                    % like 'waiting' or 'aa_doprocessing'
                    obj.workerstatus{workerIndex}=workerData{2};
                case 'done_aa_doprocessing'
                    % workerData{2} is the job index
                    jobind=workerData{2};
                    obj.isJobNotRun(jobind)=false;
                    aas_log(obj.aap,false,sprintf('PARALLEL (matlab_pct): Completed %s',obj.jobqueue(jobind).doneflag));
                    % workerData{3} is streamcache_toreturn - append to stream cache
                    processstreamcache=workerData{3};
                    if ~isempty(processstreamcache)
                        if ~isfield(obj.aap.internal,'streamcache')
                            obj.aap.internal.streamcache=processstreamcache;
                        else
                            obj.aap.internal.streamcache=[obj.aap.internal.streamcache processstreamcache];
                        end
                    end
                    % Remove this as a dependency from all of the
                    % stages dependent on the stage that has just
                    % completed
                    for depoflist=1:length(obj.jobDepOf{jobind})
                        deponmask=obj.jobDepOf{jobind}(depoflist);
                        obj.jobDepOn{deponmask}(obj.jobDepOn{deponmask}==jobind)=[];
                        obj.numJobDepOn(deponmask)=obj.numJobDepOn(deponmask)-1;
                        obj.isJobReadyToGo=[obj.isJobReadyToGo deponmask(obj.numJobDepOn(deponmask)==0)];
                    end
                otherwise
                    % trigger an error for any other input, as this may
                    % signal a bug
                    aas_log(obj.aap, true, sprintf('Unrecognized message from worker %s\n',firstOutputArg)); 
            end
        end
        
        
        function [obj]=runall_spmd_master(obj,waitforalljobs)
            % Master - sends jobs that are ready to workers
            obj.receiveWorkerMessages();
            itnum=0;
            fprintf('Now trying job execution\n');
            workerwaits=[];
            
            % Main job execution list
            while any(obj.isJobNotRun)
                submitjobstart=tic;
                itnum=itnum+1;
                
                jobreadystart=tic;
                jobreadylist=obj.isJobReadyToGo;
                obj.isJobReadyToGo=[];
                obj.benchmark.jobreadystart=obj.benchmark.jobreadystart+toc(jobreadystart);
                
                for jobreadyind=length(jobreadylist):-1:1
                    jobind=jobreadylist(jobreadyind);
                    
                    % Job has no dependencies - good to go...
                    obj.jobcount=obj.jobcount+1;
                    job=obj.jobqueue(jobind);
                    
                    % Find out who can run it, wait if no-one free
                    waittime=tic;
                    while(1)
                        waitingWorkerIndex=find(strcmp('waiting',obj.workerstatus));
                        if ~isempty(waitingWorkerIndex)
                            workerwaits(end+1)=0;
                            waitingWorkerIndex=waitingWorkerIndex(1);
                            obj.workerstatus{waitingWorkerIndex}='pending';
                            break;
                        else
                            workerwaits(end+1)=1;
                        end
                        obj.receiveWorkerMessages();
                    end
                    obj.benchmark.waitForWorkers=obj.benchmark.waitForWorkers+toc(waittime);
                    
                    % Run the job
                    aas_log(obj.aap,false,sprintf('PARALLEL (matlab_pct): Submitting to worker %d job %s',waitingWorkerIndex,job.doneflag));
                    labsendstart=tic;
                    labSend({'aa_doprocessing',obj,job,jobind},waitingWorkerIndex);
                    obj.benchmark.labsend=obj.benchmark.labsend+toc(labsendstart);
                    
                    % Flag as submitted
                    obj.isJobNotRun(jobind)=false;
                end
                obj.benchmark.submitJobs=obj.benchmark.submitJobs+toc(submitjobstart);
                
                if ~waitforalljobs
                    break;
                end
                obj.receiveWorkerMessages();
            end
            
            % Close all workers as they become free
            fprintf('Closing\n');
            while any(~strcmp('closed',obj.workerstatus(2:end)))
                waitingWorkerIndex=find(strcmp('waitforrealtime',obj.workerstatus) | strcmp('idle',obj.workerstatus) | strcmp('waiting',obj.workerstatus));
                for wwind=1:length(waitingWorkerIndex)
                    fprintf('Sending close to %d\n',waitingWorkerIndex(wwind));
                    labSend({'close'},waitingWorkerIndex(wwind));
                end
                obj.receiveWorkerMessages();
                pause(0.5);
                for w=1:length(obj.workerstatus)
                    fprintf('worker %d status %s; ',w,obj.workerstatus{w});
                end
                fprintf('\n');
            end
            
            fprintf('Proportion of jobs that had to wait for worker %f\n',mean(workerwaits));
            obj.benchmark
            
            fprintf('All done\n');
            if waitforalljobs == 1
                obj.emptyqueue;
            end
            
            if obj.fatalerrors
                aas_log(obj.aap, true, 'PARALLEL (matlab_pct): Fatal errors executing jobs.','Errors');
            end
        end
        
        
        function [obj]=runall_spmd_worker(obj)
            % Enter a loop probing whether other workers are ready to send
            % data ('hang around') and, once that is the case, either
            % spring into action by running aa_doprocessing_onetask, or
            % close up, depending on the message transmitted ('be told what
            % to do').
            % Code will be run by non-master workers only.
            alldone=false;
            while ~alldone
                % report 'waiting' status
                labSend({'status','waiting'},1);
                % Wait until *any* worker is ready to send data
                while(~labProbe)
                    pause(0.01);
                end                
                % receive the data
                workerData=labReceive;
                % the first element is the string determining the action to
                % be taken
                switch(workerData{1})
                    case 'aa_doprocessing'
                        % report new status
                        labSend({'status','aa_doprocessing'},1);
                        % copy values of received data to local variables
                        obj=workerData{2};
                        job=workerData{3};
                        jobind=workerData{4};
                        aap=obj.aap;
                        % avoid any potential confusion with global var
                        % aaworker by appending _copy to variable name
                        aaworker_copy = obj.aaworker;
                        aaworker_copy.logname = fullfile(obj.matlab_pct_path,sprintf('Job%05d_diary.txt',jobind));
                        if ~isfield(aap.internal,'streamcache')
                            nstreamcache=0;
                        else
                            nstreamcache=length(aap.internal.streamcache);
                        end
                        % *here is the beef*
                        aap=aa_doprocessing_onetask(aap,job.task,job.k,job.indices,aaworker_copy);
                        % Only return new streams...
                        if isfield(aap.internal,'streamcache') && length(aap.internal.streamcache)>nstreamcache
                            streamcache_toreturn=aap.internal.streamcache(nstreamcache+1:end);
                        else
                            streamcache_toreturn=[];
                        end
                        % report new status
                        labSend({'done_aa_doprocessing',jobind,streamcache_toreturn},1);
                    case 'close'
                        labSend({'status','closed'},1);
                        alldone=true;
                    otherwise
                        % this case is only entered if a non-master worker
                        % is ready to send data, which does not affect
                        % processing of the data, but may be worth
                        % reporting (if only for debugging purposes)
                        aas_log(obj.aap, false, sprintf('Message from worker: %s\n',workerData{1}));                         
                end
            end
        end
        
        
        function [obj]=waitforrealtime(obj)
            % THREAD TO ROUTE INCOMING REALTIME DATA
            % Overview:
            %  A set of realtime modules (like aamod_realtime_wait_epi)
            %  contain the attributes waitforrealtime_singlefile or
            %  (for future implementation) waitforrealtime_wholeseries
            %  As these jobs are added to the matlab_pct queue, they
            %  are given (negative) dependencies that will never be satistifed by the
            %  normal aa pathway.
            % When realtime metadata arrives in the series*.txt files, then
            %  the following code attempts to match it with the waiting
            %  realtime modules. When it finds a match, it creates the
            %  output and marks the module as complete.
            %
            % Incoming realtime data
            
            impth=fullfile(obj.aap.directory_conventions.realtime.path,'incomingmeta');
            imfid=[];
            filenumber=0;
            lastdatatime=tic;
            status='waitforrealtime';
            labSend({'status','waitforrealtime'},1);
            while(1)
                if labProbe
                    [data source tag]=labReceive;
                    switch(data{1})
                        case 'close'
                            labSend({'status','closed'},1);
                            break;
                    end
                end
                if isempty(imfid)
                    allseries=dir(fullfile(impth,'series*.txt'));
                    while isempty(imfid)
                        if filenumber<length(allseries)
                            %
                            % Open next series file
                            filenumber=filenumber+1;
                            fprintf('Going to look for metafile %s\n',allseries(filenumber).name);
                            imfid=fopen(fullfile(impth,allseries(filenumber).name),'r');
                            linesread=0;
                            % Trap for empty series files
                            if feof(imfid)
                                fclose(imfid);
                                imfid=[];
                            end
                            if strcmp(status,'idle')
                                labSend({'status','waitforrealtime'},1);
                                status='waitforrealtime';
                            end
                        else
                            % More than 5s since the last data?
                            if toc(lastdatatime)>5
                                labSend({'status','idle'},1);
                                status='idle';
                                pause(0.5);
                            end
                            break;
                        end
                    end
                end
                
                % Read a line from it and process
                if ~isempty(imfid)
                    % Reopen and throw away previously read lines
                    imfid=fopen(fullfile(impth,allseries(filenumber).name),'r');
                    %                     fprintf('Getting ready to junk %d lines\n',linesread);
                    for junklines=1:linesread
                        ln=fgetl(imfid);
                    end
                    %                     fprintf('Done junking\n');
                    
                    % Any new lines?
                    if ~feof(imfid)
                        lastdatatime=tic;
                        ln=strtrim(fgetl(imfid));
                        linesread=linesread+1;
                        [cmd rem]=strtok(ln,[' ' 9 10]);
                        cmd=strtrim(cmd);
                        fprintf('got line %s\n',ln);
                        if strcmp(cmd,'MEAS_START')
                            aas_log(obj.aap,false,sprintf('Realtime incoming %s: MEAS_START',allseries(filenumber).name));
                            if ~feof(imfid)
                                ln=strtrim(fgetl(imfid));
                                linesread=linesread+1;
                                [cmd rem]=strtok(ln,[' ' 9 10]);
                                cmd=strtrim(cmd);
                            end
                        end
                        
                        if strcmp(cmd,'MEAS_FINISHED')
                            aas_log(obj.aap,false,sprintf('Realtime incoming %s: MEAS_FINISHED',allseries(filenumber).name));
                            fclose(imfid);
                            imfid=[];
                        elseif strcmp(cmd,'DATAFILE')
                            % May need to enclose the following in a
                            % try-catch in case file is only partly
                            % written
                            [datatype rem]=strtok(rem,[' ' 9 10]);
                            datatype=strtrim(datatype);
                            switch (datatype)
                                case 'DICOMRAW'
                                    
                                case 'DICOMIMA'
                                    
                                    if linesread==3 || linesread>23
                                        dcmfile=fullfile(obj.aap.directory_conventions.realtime.path,'rawdata',strtrim(rem))
                                        dcmfile(dcmfile=='\')='/';
                                        
                                        aas_log(obj.aap,false,sprintf('Realtime incoming %s: ',allseries(filenumber).name,ln));
                                        
                                        % Throw a strop, for now at least
                                        if ~exist(dcmfile,'file')
                                            aas_log(obj.aap,true,sprintf('Dicom file not found: %s',dcmfile))
                                        end
                                        
                                        % Load up DICOM and see what it will
                                        % satisfy
                                        H=spm_dicom_headers(dcmfile);
                                        
                                        fprintf('ImageType is %s and ProtocolName %s\n',H{1}.ImageType,H{1}.ProtocolName)
                                        % Find out if we can satisfy any waiting
                                        % dependencies with this
                                        if isfield(obj.realtime_deps,'satisfied') && ~isempty(obj.realtime_deps)
                                            rdlist=find(~[obj.realtime_deps.satisfied]);
                                            for rdind=1:length(rdlist)
                                                rdtarget=obj.realtime_deps(rdlist(rdind));
                                                fprintf('Looking for something where %s is %s\n',rdtarget.dcmfield,strtrim(H{1}.(rdtarget.dcmfield)));
                                                if strcmp(strtrim(H{1}.(rdtarget.dcmfield)),strtrim(rdtarget.dcmfilter))
                                                    % Bingo, found the job we're going
                                                    % to satisfy.
                                                    tic
                                                    dcmjob=obj.jobqueue(rdtarget.njob);
                                                    obj.aap=aas_setcurrenttask(obj.aap,dcmjob.k);
                                                    % Fake its output data
                                                    domainpth=aas_getpath_bydomain(obj.aap,obj.jobqueue(rdtarget.njob).domain,obj.jobqueue(rdtarget.njob).indices);
                                                    aas_makedir(obj.aap,domainpth);
                                                    [pth nme ext]=fileparts(dcmfile);
                                                    dcmfile_jobpth=fullfile(domainpth,[nme ext]);
                                                    copyfile(dcmfile,dcmfile_jobpth);
                                                    
                                                    % Register the output stream
                                                    obj.aap=aas_desc_outputs(obj.aap,dcmjob.domain,dcmjob.indices,obj.aap.schema.tasksettings.(dcmjob.stagename)(obj.aap.tasklist.main.module(dcmjob.k).index).outputstreams.stream{1}, dcmfile_jobpth);
                                                    
                                                    % This dependency is done
                                                    obj.realtime_deps(rdlist(rdind)).satisfied=true;
                                                    
                                                    obj.isJobNotRun(rdtarget.njob)=false;
                                                    
                                                    % Write done flag
                                                    doneflag=aas_doneflag_getpath_bydomain(obj.aap,dcmjob.domain,dcmjob.indices,dcmjob.k);
                                                    aas_writedoneflag(obj.aap,doneflag);
                                                    
                                                    % Send a message to say this has
                                                    % been completed, providing only
                                                    % the last streamcache entry
                                                    labSend({'done_aa_doprocessing',rdtarget.njob,obj.aap.internal.streamcache(end)},1);
                                                    
                                                    % FOR SIMULATIONS ONLY: Max rate of new data
                                                    %                                                 fprintf('!!!!!!!!!!!!! Simulation pause\n');
                                                    %                                                 pause(0.5);
                                                    break;
                                                else
                                                    fprintf('It is %s and doesn''t match %s\n',strtrim(H{1}.(rdtarget.dcmfield)),strtrim(rdtarget.dcmfilter));
                                                end
                                            end
                                        end
                                    end
                                otherwise
                                    aas_log(obj.aap,true,sprintf('Unknown realtime incoming datatype %s in %s\n',datatype,allseries(filenumber).name));
                            end
                        end
                    end
                    if ~isempty(imfid)
                        fclose(imfid);
                    end
                
                end
                
                
            end
        end
    end
end