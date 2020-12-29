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
    % One worker is defined as the master; communication between master and
    % the workers is accomplished via labSend, labReceive, and labProbe.
    %
    % The sequence of method/function calls triggering the execution of a
    % module's code is
    %
    %   runall > runall_spmd_master | runall_spmd_worker > aa_doprocessing_onetask
    % 
    %
    % aaq_matlab_pct Properties (*inherited):
    %   toHandlePool - logical, false indicating that superclass handles
    %                  parpool
    %   benchmark    - struct, collects job statistics
    %   matlab_pct_path - char arr, subdirectory for storing job diaries
    %   jobcount     - double, counter for jobs in pipeline
    %   *aap         - struct
    %   *isOpen      - logical, flag indicating whether Taskqueue is open
    %   *fatalerrors - logical, flag indicating fatal error
    %   *jobqueue    - struct array, composed of 'taskmasks'
    % 
    % aaq_matlab_pct Methods (*inherited):
    %   addtask             - Call superclass' addtask and additionally scan
    %                         for dependencies
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
    %   *save              - save self to file
    %   *emptyqueue        - clear the task queue (set jobqueue to [])
    %   *allocate          - allocate a job from the task queue to a worker
    %   *getjobdescription - return struct task

    
    properties
        toHandlePool    = false; % aaq handles pools
        benchmark       = struct('receiveWorkerMessages',0,...
                            'submitJobs',0,...
                            'waitForWorkers',0,...
                            'labsend',0,...
                            'jobreadystart',0);
        matlab_pct_path = []; % subdirectory for storing job diaries (<aaworker.parmpath>/matlab_pct)
        jobcount        = 0;  % counter for jobs in pipeline
    end
    
    properties (Hidden)
        jobstudypths            = {} % actualised studypath (for each job)
        readytogo               = [] % green sign (for each job)
        jobnotrun               = [] % job status (for each job)

        dep_names               = {} % doneflag (for each job)
        dep_done                = [] % is doneflag exists (for each job)
        depon                   = {} % jobs the actual job depends on (for each job)
        depof                   = {} % jobs depending on the actual job (for each job)
        depon_num               = [] % # jobs the actual job depends on (for each job)

        % real-time
        realtime_deps           = []

        % low-level
        aaworker                = [] % aaworker struct
        aaparallel              = [] % aaparallel struct
        workerstatus            = {} % worker status (for each worker)
        
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
                obj.toHandlePool = true;                
            end
        end

        
        function close(obj)
            % Closes matlabpool and invokes super class' close function
            if obj.toHandlePool
                aas_matlabpool('close'); 
            end
            close@aaq(obj);
        end
        
        
        function obj=addtask(obj,taskmask)
            % Add a task to the task queue.
            % Adds to the parent class's method by scanning the
            % dependencies of the new module.
            obj=addtask@aaq(obj,taskmask);
            index=obj.aap.tasklist.main.module(taskmask.k).index;
            njob=length(obj.jobqueue);
            obj.jobnotrun(njob)=true;
            
            % For better performance, make a list of all dependencies and
            % refer to them by index 
            % Removes need for endless done flag file checking
            obj.depof{njob}=[];
            obj.depon{njob}=[];
            obj.depon_num(njob)=0;
            
            % Store this task as a future dependency of others
            obj.dep_names{njob}=taskmask.doneflag;
            obj.dep_done(njob)=exist(taskmask.doneflag,'file');
            obj.jobstudypths{njob}=aas_getstudypath(obj.aap,taskmask.k);
            
            % Does this stage need realtime input?
            if isfield(obj.aap.schema.tasksettings.(taskmask.stagename)(index).ATTRIBUTE,'waitforrealtime_singlefile') && ~isempty(obj.aap.schema.tasksettings.(taskmask.stagename)(index).ATTRIBUTE.waitforrealtime_singlefile)
                [dcmfield, dcmfilter]=strtok(obj.aap.schema.tasksettings.(taskmask.stagename)(index).ATTRIBUTE.waitforrealtime_singlefile,'=');
                newrtd=struct('dcmfield',dcmfield,'dcmfilter',strtrim(dcmfilter(2:end)),'njob',njob,'eventtype','singlefile','satisfied',false);
                if isempty(obj.realtime_deps)
                    obj.realtime_deps=newrtd;
                else
                    obj.realtime_deps(end+1)=newrtd;
                end
                obj.depon{njob}=-length(obj.realtime_deps); % minus indicates realtime dependency
                obj.depon_num(njob)=1;
            end
            
            % Go through each dependency of this task and find the index of
            % that
            for tbcfind=1:length(obj.jobqueue(njob).tobecompletedfirst)
                % Find or put this dependency in a list
                df=obj.jobqueue(njob).tobecompletedfirst{tbcfind};
                depind=find(strcmp(df,obj.dep_names));
                if isempty(depind)
                    aas_log(aap,true,sprintf('Dependency %s of %s not found.',df,taskmask.doneflag));
                end
                
                % Only record dependencies that haven't already been
                % satisfied
                if ~obj.dep_done(depind)
                    obj.depon{njob}(end+1)=depind;
                    obj.depon_num(njob)=obj.depon_num(njob)+1;
                    obj.depof{depind}(end+1)=njob;
                end
            end
            
            % Add it to the readytogo queue if appropriate
            if obj.depon_num(njob)==0
                obj.readytogo=[obj.readytogo njob];
            end
        end

        
        function [obj]=runall(obj, ~, waitforalljobs)
            % Run all jobs, using matlabs "single program, multiple data"
            % construct. This isn't really what we use it for - we actually
            % set pass messages between the modules to make a queue happen.
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
            % Used by the master thread to pump any messages from workers
            tic
            messagesreceived=false;
            while(labProbe)
                [data, source, ~]=labReceive;
                obj.processWorkerMessage(data,source);
                messagesreceived=true;
            end
            obj.benchmark.receiveWorkerMessages=obj.benchmark.receiveWorkerMessages+toc;
        end
        
        
        function [obj]=processWorkerMessage(obj,data,source)
            % Execute this worker message
            switch(data{1})
                case 'status'
                    obj.workerstatus{source}=data{2};
                case 'done_aa_doprocessing'
                    jobind=data{2};
                    obj.jobnotrun(jobind)=false;
                    aas_log(obj.aap,false,sprintf('PARALLEL (matlab_pct): Completed %s',obj.jobqueue(jobind).doneflag));
                    
                    % Append to stream cache with results from this process
                    processstreamcache=data{3};
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
                    for depoflist=1:length(obj.depof{jobind})
                        deponmask=obj.depof{jobind}(depoflist);
                        obj.depon{deponmask}(obj.depon{deponmask}==jobind)=[];
                        obj.depon_num(deponmask)=obj.depon_num(deponmask)-1;
                        obj.readytogo=[obj.readytogo deponmask(obj.depon_num(deponmask)==0)];
                    end
                otherwise
                    fprintf('Unrecognized message from worker %s\n',data{1});
                    
            end
        end
        
        
        function [obj]=runall_spmd_master(obj,waitforalljobs)
            % Master - sends jobs that are ready to workers
            obj.receiveWorkerMessages();
            itnum=0;
            fprintf('Now trying job execution\n');
            workerwaits=[];
            
            % Main job execution list
            while any(obj.jobnotrun)
                submitjobstart=tic;
                itnum=itnum+1;
                
                jobreadystart=tic;
                jobreadylist=obj.readytogo;
                obj.readytogo=[];
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
                    obj.jobnotrun(jobind)=false;
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
            
            if (obj.fatalerrors)
                aas_log(obj.aap,false,'PARALLEL (matlab_pct): Fatal errors executing jobs. The errors were:','Errors');
                for errind=1:length(errline)
                    aas_log(obj.aap,false,[' ' errline{errind}],'Errors');
                end
                aas_log(obj.aap,true,'PARALLEL (matlab_pct): Terminating because of errors','Errors');
            end
        end
        
        
        function [obj]=runall_spmd_worker(obj)
            % Code for a worker
            % Basically hangs around waiting to be told what to do, and
            % then executes code when asked
            alldone=false;
            while ~alldone
                % Say we're waiting for something to do
                labSend({'status','waiting'},1);
                
                % Wait to be told what that is
                while(~labProbe)
                    pause(0.01);
                end                
                
                [data, ~, ~]=labReceive;
                switch(data{1})
                    case 'aa_doprocessing'
                        labSend({'status','aa_doprocessing'},1);
                        obj=data{2};
                        job=data{3};
                        jobind=data{4};
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
                        aap=aa_doprocessing_onetask(aap,job.task,job.k,job.indices,aaworker_copy);
                        % Only return new streams...
                        if isfield(aap.internal,'streamcache') && length(aap.internal.streamcache)>nstreamcache
                            streamcache_toreturn=aap.internal.streamcache(nstreamcache+1:end);
                        else
                            streamcache_toreturn=[];
                        end
                        labSend({'done_aa_doprocessing',jobind,streamcache_toreturn},1);
                    case 'close'
                        labSend({'status','closed'},1);
                        alldone=true;
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
                                                    
                                                    obj.jobnotrun(rdtarget.njob)=false;
                                                    
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