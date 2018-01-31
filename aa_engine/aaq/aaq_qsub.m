classdef aaq_qsub<aaq
    properties
        pool = []
        QV = []
    end
    properties (Hidden)
        jobnotrun = []
        jobinfo = []
        jobretries = []
        waitforalljobs
        
        % ensure MAXFILTER license
        initialSubmitArguments = '';
    end
    methods
        function [obj]=aaq_qsub(aap)
            global aaworker;
            global aaparallel;
            
            try
                if ~isempty(aap.directory_conventions.poolprofile)
                    % Parse configuration
                    queue = '';
                    poolprofile = textscan(aap.directory_conventions.poolprofile,'%s','delimiter', ':'); poolprofile = poolprofile{1};
                    if numel(poolprofile) == 1, poolprofile = poolprofile{1};
                    else [poolprofile, queue] = poolprofile{1:2}; end
                    
                    profiles = parallel.clusterProfiles;
                    if ~any(strcmp(profiles,poolprofile))
                        ppfname = which(spm_file(poolprofile,'ext','.settings'));
                        if isempty(ppfname)
                            aas_log(obj.aap,true,sprintf('ERROR: settings for pool profile %s not found!',poolprofile));
                        else
                            obj.pool=parcluster(parallel.importProfile(ppfname));
                        end
                    else
                        aas_log(obj.aap,false,sprintf('INFO: pool profile %s found',poolprofile));
                        obj.pool=parcluster(poolprofile);
                    end
                    switch class(obj.pool)
                        case 'parallel.cluster.Torque'
                            aas_log(obj.aap,false,'INFO: pool Torque is detected');
                            obj.pool.ResourceTemplate = sprintf('-l nodes=^N^,mem=%dGB,walltime=%d:00:00', aaparallel.memory,aaparallel.walltime);
                            if any(strcmp({aap.tasklist.main.module.name},'aamod_meg_maxfilt')) && ... % maxfilt module detected
                                    ~isempty(aap.directory_conventions.neuromagdir) % neuromag specified
                                obj.initialSubmitArguments = ' -W x=\"NODESET:ONEOF:FEATURES:MAXFILTER\"';
                            end
                            obj.pool.SubmitArguments = strcat(obj.pool.SubmitArguments,obj.initialSubmitArguments);
                            aaparallel.numberofworkers = 1;
                        case 'parallel.cluster.LSF'
                            aas_log(obj.aap,false,'INFO: pool LSF is detected');
                            if ~isempty(queue), obj.initialSubmitArguments = [' -q ' queue]; end
                            obj.initialSubmitArguments = sprintf('%s -M %d -R "rusage[mem=%d]"',obj.initialSubmitArguments,aaparallel.memory*1000,aaparallel.memory*1000);
                            obj.pool.SubmitArguments = strcat(obj.pool.SubmitArguments,obj.initialSubmitArguments);
                            aaparallel.numberofworkers = aap.options.aaparallel.numberofworkers;
                        case 'parallel.cluster.Generic'
                            aas_log(obj.aap,false,'INFO: Generic engine is detected');
                            obj.pool.IndependentSubmitFcn = obj.SetArg(obj.pool.IndependentSubmitFcn,'walltime',aaparallel.walltime);
                            obj.pool.IndependentSubmitFcn = obj.SetArg(obj.pool.IndependentSubmitFcn,'memory',aaparallel.memory);
                            aaparallel.numberofworkers = 1;
                        case 'Local'
                            aas_log(obj.aap,false,'INFO: Local engine is detected');  
                            aaparallel.numberofworkers = aap.options.aaparallel.numberofworkers;
                    end
                else
                    obj.pool = parcluster('local');
                    aaparallel.numberofworkers=12;
                end
                obj.pool.NumWorkers = aaparallel.numberofworkers;
                obj.pool.JobStorageLocation = aaworker.parmpath;
            catch ME
                aas_log(aap,false,'WARNING: Cluster computing is not supported!');
                aas_log(aap,false,sprintf('\tERROR in %s:\n\tline %d: %s',ME.stack(1).file, ME.stack(1).line, ME.message),aap.gui_controls.colours.warning);
                obj.pool=[];
            end
            obj.aap=aap;
        end
        
        function close(obj)
            if ~isempty(obj.pool)
                for j = 1:numel(obj.pool.Jobs)
                    obj.pool.Jobs(j).cancel;
                end
            end
            close@aaq(obj);
        end
        
        %% Queue jobs on Qsub:
        %  Queue job
        %  Watch output files
        
        % Run all tasks on the queue
        function [obj]=runall(obj, dontcloseexistingworkers, waitforalljobs) %#ok<INUSL>
            obj.waitforalljobs = waitforalljobs;
            
            global aaworker
            
            % Check number of jobs & monitored files
            njobs=length(obj.jobqueue);
            
            % We have already submitted some of these jobs
            submittedJobs = 1:length(obj.jobnotrun);
            obj.jobnotrun = true(njobs,1);
            obj.jobnotrun(submittedJobs) = false;
                
            % Create a array of job retry attempts (needs to be specific to
            % each module so using dynamic fields with modulename as fieldname. 
            % modulename is also saved in jobinfo so is easily retrieved later)
            if ~isempty(obj.jobqueue)
                % probably a better place to get current module index.
                modulename = obj.aap.tasklist.main.module(obj.jobqueue(end).k).name; 
                obj.jobretries.(modulename) = ones(njobs,1);
                
                if obj.jobqueue(end).k == length(obj.aap.tasklist.main.module)
                    obj.waitforalljobs = true;
                end
            end
            
            jobqueuelimit = obj.aap.options.aaparallel.numberofworkers;
            printswitches.jobsinq = true; % switches for turning on and off messages
            
            while any(obj.jobnotrun) || (obj.waitforalljobs && ~isempty(obj.jobinfo))
                % Lets not overload the filesystem
                pause(0.1);
                
                pool_length = length(obj.pool.Jobs);
                nfreeworkers = jobqueuelimit - pool_length;
                
                % Only run if there are free workers, otherwise display
                % message that no workers are available
                if nfreeworkers > 0
                    printswitches.nofreeworkers = true; % reset the display command
                    % Find how many free workers available, then allocate next
                    % batch. Skip section if there are no jobs to run.
                    if any(obj.jobnotrun)
                        if printswitches.jobsinq;
                            aas_log(obj.aap, false, sprintf('Jobs in aa queue: %d\n', sum(obj.jobnotrun)))
                            printswitches.jobsinq = false; % Don't display again unless queue length changes
                        end
                        runjobs = shiftdim(find(obj.jobnotrun))';
                        nfreeworkers = min([nfreeworkers, length(runjobs)]);
                        runjobs = runjobs(1:nfreeworkers);
                        readytorunall = true(size(runjobs));
                        for i = runjobs
                            rtrind = runjobs == i; % index for readytorunall
                            if (obj.jobnotrun(i))
                                % Find out whether this job is ready to be allocated by
                                % checking dependencies (done_ flags)
                                for j=1:length(obj.jobqueue(i).tobecompletedfirst)
                                    if (~exist(obj.jobqueue(i).tobecompletedfirst{j},'file'))
                                        readytorunall(rtrind) = false;
                                    end
                                end
                                
                                if readytorunall(rtrind)
                                    % Add the job to the queue, and create
                                    % the job info in obj.jobinfo
                                    obj.add_from_jobqueue(i);
                                    printswitches.jobsinq = true;
                                end
                            end
                        end
                        
                        if ~any(readytorunall)
                            % display monitor and update job states
                            aas_log(obj.aap, false, 'Workers available, but no jobs ready to run. Waiting 60 seconds...')
                            obj.job_monitor(true);
                            pause(60)
                        else
                            % silently update job states
                            obj.job_monitor(true);
                        end
                    elseif ~isempty(obj.jobinfo)
                        % If no jobs left then monitor and update job states
                        aas_log(obj.aap, false, 'No jobs in the queue, waiting for remaining jobs to complete...')
                        obj.job_monitor_loop;
                        pause(60)
                    end
                else
                    aas_log(obj.aap, false, 'No free workers available: monitoring the queue...')
                    obj.job_monitor_loop;
                end
                
                idlist = [obj.jobinfo.JobID];
                for id = idlist
                    JI = obj.jobinfo([obj.jobinfo.JobID] == id); 
                    % For clarity use JI.JobID from now on (JI.JobID = id). All
                    % job information is stored in JI, including the main
                    % queue index JI.qi used to refer back to the original
                    % job queue (obj.jobqueue) created when this object was called.
                    
                    % This section of code only runs if the job is not
                    % running (has finished or error)
                    if strcmp(JI.state, 'finished') || strcmp(JI.state, 'error') || strcmp(JI.state, 'failed')
                        Jobs = obj.pool.Jobs([obj.pool.Jobs.ID] == JI.JobID);
                        if isempty(Jobs) % cleared by the GUI
                            if obj.QV.isvalid
                                obj.QV.Hold = false;
                            end
                            obj.fatalerrors = true; % abnormal terminations
                            obj.close;
                            return;
                        end
                        
                        %                     moduleName = JI.InputArguments{1}.tasklist.main.module(JI.InputArguments{3}).name;
                        aap = aas_setcurrenttask(obj.aap,JI.InputArguments{3});
                        moduleName = obj.aap.tasklist.main.module(JI.InputArguments{3}).name;
                        indices = JI.InputArguments{4};
                        datname = ''; datpath = '';
                        if numel(indices) > 0 % subject specified
                            datname = aas_getsubjdesc(aap,indices(1));
                            datpath = aas_getsubjpath(aap,indices(1));
                        end
                        if numel(indices) > 1 % session specified
                            datname = aas_getsessdesc(aap,indices(1),indices(2));
                            datpath = aas_getsesspath(aap,indices(1),indices(2));
                        end
                                                
                        if ~isempty(Jobs.Tasks.Error)
                            switch Jobs.Tasks.Error.identifier
                                case 'parallel:job:UserCancellation'
                                    state = 'cancelled';
                                otherwise
                                    state = 'error';
                            end
                        end
                        
                        switch JI.state
                            case 'pending'
                                t = toc(JI.tic);
                                % aa to switch this on/off or extend the time? On very busy
                                % servers this might cause all jobs to be
                                % perpetually deleted and restarted.
                                if (obj.aap.options.aaworkermaximumretry > 1) && (t > obj.aap.options.aaworkerwaitbeforeretry) % if job has been pending for more than N seconds
                                    obj.remove_from_jobqueue(JI.JobID, true); % 2nd argument = retry
                                end
                                
                            case 'failed' % failed to launch
                                msg = sprintf(...
                                    ['Failed to launch (Licence?)!\n'...
                                    'Check <a href="matlab: open(''%s'')">logfile</a>\n'...
                                    'Queue ID: %d | qsub ID %d | Subject ID: %s'],...
                                    fullfile(obj.pool.JobStorageLocation, Jobs.Tasks.Parent.Name, [Jobs.Tasks.Name '.log']),...
                                    JI.qi, JI.JobID, JI.subjectinfo.subjname);
                                % If there is an error, it is fatal...
                                
                                fatalerror = obj.jobretries.(JI.modulename)(JI.qi) > obj.aap.options.aaworkermaximumretry; % will be true if max retries reached
                                if fatalerror
                                    msg = sprintf('%s\nMaximum retries reached\n', msg);
                                    aas_log(obj.aap,fatalerror,msg,obj.aap.gui_controls.colours.error)
                                end
                                
                                % This won't happen if the max retries is
                                % reached, allowing debugging
                                obj.remove_from_jobqueue(JI.JobID, true); % 2nd argument = retry
                                    
                            case 'cancelled' % cancelled
                                msg = sprintf('Job%d had been cancelled by user!\n Check <a href="matlab: open(''%s'')">logfile</a>\n',JI.JobID,...
                                    fullfile(obj.pool.JobStorageLocation,Jobs.Tasks.Parent.Name,[Jobs.Tasks.Name '.log']));
                                % If there is an error, it is fatal...
                                aas_log(obj.aap,true,msg,obj.aap.gui_controls.colours.warning)
                                
                            case 'finished' % without error
                                if isempty(Jobs.Tasks.FinishTime), continue; end
                                msg = sprintf('JOB %d: \tMODULE %s \tON %s \tSTARTED %s \tFINISHED %s \tUSED %s.',...
                                    JI.JobID,moduleName,datname,Jobs.Tasks.CreateTime,Jobs.Tasks.FinishTime,...
                                    sprintf('%*dd %*dh %dm %ds',round(datevec(jts_etime(Jobs.Tasks.FinishTime,Jobs.Tasks.CreateTime)/3600/24))));
                                aas_log(obj.aap,false,msg,obj.aap.gui_controls.colours.completed);
                                
                                % Also save to file with module name attached!
                                fid = fopen(fullfile(aaworker.parmpath,'qsub','time_estimates.txt'), 'a');
                                fprintf(fid,'%s\n',msg);
                                fclose(fid);
                                
                                obj.remove_from_jobqueue(JI.JobID, false);
                                
                            case 'error' % running error
                                
                                % Check whether the error was a "file does not
                                % exist" type. This can happen when a dependent
                                % folder is only partially written upon job execution.
                                % jobinfo etc is indexed by the job ID so get i from jobinfo.
                                
                                if obj.jobretries.(JI.modulename)(JI.qi) <= obj.aap.options.aaworkermaximumretry
                                    msg = sprintf(['%s\n\n JOB FAILED WITH ERROR: \n %s',...
                                        ' \n\n Waiting 60 seconds then trying again',...
                                        ' (%d tries remaining for this job)\n'...
                                        'Press Ctrl+C now to quit, then run aaq_qsub_debug()',...
                                        ' to run the job locally in debug mode.\n'],...
                                        Jobs.Tasks.Diary, Jobs.Tasks.ErrorMessage, obj.aap.options.aaworkermaximumretry - obj.jobretries.(JI.modulename)(JI.qi));
                                    aas_log(aap, false, msg);
                                    obj.jobnotrun(JI.qi) = true;
                                    obj.remove_from_jobqueue(JI.JobID, true);
                                    pause(60)
                                else
                                    msg = sprintf('Job%d on <a href="matlab: cd(''%s'')">%s</a> had an error: %s\n',JI.JobID,datpath,datname,Jobs.Tasks.ErrorMessage);
                                    for e = 1:numel(Jobs.Tasks.Error.stack)
                                        % Stop tracking to internal
                                        if strfind(Jobs.Tasks.Error.stack(e).file,'distcomp'), break, end
                                        msg = [msg sprintf('<a href="matlab: opentoline(''%s'',%d)">in %s (line %d)</a>\n', ...
                                            Jobs.Tasks.Error.stack(e).file, Jobs.Tasks.Error.stack(e).line,...
                                            Jobs.Tasks.Error.stack(e).file, Jobs.Tasks.Error.stack(e).line)];
                                    end
                                    % If there is an error, it is fatal...
                                    obj.fatalerrors = true;
                                    aas_log(obj.aap,true,msg,obj.aap.gui_controls.colours.error)
                                    obj.jobnotrun(JI.qi) = true;
                                end
                        end
                    end
                end
               
                obj.QVUpdate; 
            end

        end
        
        function obj = QVUpdate(obj)
            if obj.aap.options.aaworkerGUI
                % queue viewer
                if ~isempty(obj.pool)
                    if ~isempty(obj.QV) && ~obj.QV.isvalid % started but killed
                        return
                    end
                    if (isempty(obj.QV) || ~obj.QV.OnScreen) % not started or closed
                        obj.QV = aas_qsubViewerClass(obj);
                        obj.QV.Hold = true;
                        obj.QV.setAutoUpdate(false);
                    else
                        obj.QV.UpdateAtRate;
                        if obj.waitforalljobs, obj.QV.Hold = false; end
                    end
                end
            end
        end
        
        function obj = QVClose(obj)
            if obj.aap.options.aaworkerGUI
                if ~isempty(obj.QV) && obj.QV.isvalid
                    obj.QV.Close;
                    obj.QV.delete;
                    obj.QV = [];
                end
            end
        end
        
        function obj = pool_args(obj,varargin)
            global aaparallel;
            memory = aaparallel.memory;
            walltime = aaparallel.walltime;
            
            for iarg = 1:numel(varargin)
                if ~ischar(varargin{iarg}), continue; end
                switch varargin{iarg}
                    case 'mem'
                        if ~isempty(varargin{iarg+1}), memory = varargin{iarg+1}; end
                    case 'walltime'
                        if ~isempty(varargin{iarg+1}), walltime = varargin{iarg+1}; end
                end
            end
            
            switch class(obj.pool)
                case 'parallel.cluster.Torque'
                    if round(memory) == memory % round
                        memory = sprintf('%dGB',memory);
                    else % non-round --> MB
                        memory = sprintf('%dMB',memory*1000);
                    end
                    obj.pool.SubmitArguments = strcat(sprintf('-q compute -l mem=%s -l walltime=%d',memory,walltime*3600),obj.initialSubmitArguments);
                    %                 obj.pool.SubmitArguments = strcat(obj.initialSubmitArguments,...
                    %                     sprintf(' -N Mod%02d_',job.k),...
                    %                     sprintf('%03d',job.indices));
                case 'parallel.cluster.Generic'
                    aas_log(aap,false,'INFO: Generic engine is detected');
                    obj.pool.IndependentSubmitFcn = obj.SetArg(obj.pool.IndependentSubmitFcn,'walltime',walltime);
                    obj.pool.IndependentSubmitFcn = obj.SetArg(obj.pool.IndependentSubmitFcn,'memory',memory);
            end
        end
        
        
        function [obj]=qsub_q_job(obj,job)
            global aaworker
            global aacache
            % Let's store all our qsub thingies in one particular directory
            qsubpath=fullfile(aaworker.parmpath,'qsub');
            aas_makedir(obj.aap,qsubpath);
            cd(qsubpath);
            % Submit the job
            if ~isempty(obj.pool)
                % Check how much memory and time we should assign to the job
                qsubsettings = {'mem',[],'walltime',[]};
                % if isfield(obj.aap.tasksettings.(job.stagename)(obj.aap.tasklist.main.module(job.k).index),'qsub')
                % qsub = obj.aap.tasksettings.(job.stagename)(obj.aap.tasklist.main.module(job.k).index).qsub;
                % for f = fieldnames(qsub)'
                % switch f{1}
                % case 'memoryBase'
                % qsubsettings{2} = qsub.memoryBase;
                % case 'timeBase'
                % qsubsettings{4} = qsub.timeBase;
                % end
                % end
                % end
                if isa(obj.pool,'parallel.cluster.Torque'), obj = obj.pool_args(qsubsettings{:}); end
                J = createJob(obj.pool);
                cj = @aa_doprocessing_onetask;
                nrtn = 0;
                inparg = {obj.aap,job.task,job.k,job.indices, aaworker};
                if isprop(J,'AutoAttachFiles'), J.AutoAttachFiles = false; end
                % [RT 2013-09-04 and 2013-11-11; TA 2013-11-14 and 2014-12-12] Make workers self-sufficient by passing
                % them the aa paths. Users don't need to remember to update
                % their own default paths (e.g. for a new aa version)
                if isprop(J,'AdditionalPaths')
                    J.AdditionalPaths = aacache.path.reqpath;
                elseif isprop(J,'PathDependencies')
                    J.PathDependencies = aacache.path.reqpath;
                end
                createTask(J,cj,nrtn,inparg,'CaptureDiary',true);
                success = false;
                retries = 0;
                % Job submission can sometimes fail (server fault) (DP). Added rety to cope this this.
                while success == false
                    try
                        J.submit;
                        success = true;
                    catch ME
                        if retries > obj.aap.options.aaworkermaximumretry
                            throw(ME)
                        else
                            aas_log(obj.aap, false, 'Could not add job. Retrying...')
                            retries = retries + 1;
                            pause(5)
                        end
                    end
                end
                % % State what the assigned number of hours and GB is...
                % Naas_movParsot in use [TA]
                % fprintf('Job %s, assigned %0.4f hours. and %0.9f GB\n\n', ...
                % job.stagename, timReq./(60*60), memReq./(1024^3))
            else
                aa_doprocessing_onetask(obj.aap,job.task,job.k,job.indices);
            end
        end
        
        
        function obj = add_from_jobqueue(obj, i)
            global aaworker
            try
                % Add a job to the queue
                job=obj.jobqueue(i);
                obj.qsub_q_job(job);

                % Create job info for referencing later
                % (also clean up done jobs to prevent IDs occuring twice)
                latestjobid = max([obj.pool.Jobs.ID]);
                if ~all(obj.jobnotrun(i)) % if any jobs have been run yet
                    obj.jobinfo([obj.jobinfo.JobID] == latestjobid) = []; % remove prev job with same ID
                end
                
                ji.InputArguments = {[], job.task, job.k, job.indices, aaworker};
                ji.modulename = obj.aap.tasklist.main.module(ji.InputArguments{3}).name;
                [junk, ji.jobpath]=aas_doneflag_getpath_bydomain(obj.aap,job.domain,job.indices,job.k);
                ji.JobID = latestjobid;
                ji.qi = i;
                ji.jobrunreported = false;
                ji.state = 'pending';
                
                if strcmp(job.domain, 'study')
                    ji.subjectinfo = struct('subjname', 'ALL SUBJECTS');
                else
                    ji.subjectinfo = obj.aap.acq_details.subjects(job.indices(1));
                end
                
                obj.jobinfo = [obj.jobinfo, ji];
                obj.jobnotrun(i) = false;
                rt = obj.jobretries.(ji.modulename)(ji.qi);
                if isnan(rt), rt = 0; end
                aas_log(obj.aap, false, sprintf('Added job %s with qsub ID %3.1d | Subject ID: %s | Retry: %3.1d | Jobs submitted: %3.1d',...
                    ji.modulename, ji.JobID, ji.subjectinfo.subjname, rt, length(obj.pool.Jobs)))
            catch ME
                if strcmp(ME.message, 'parallel:job:OperationOnlyValidWhenPending')
                    obj.jobretries.(ji.modulename)(i) = obj.jobretries.(ji.modulename)(i) + 1;
                    aas_log(obj.aap, false, sprintf('WARNING: Error starting job: %s | Retries: %d', ME.message, obj.jobretries.(ji.modulename)(i)))
                    obj.jobnotrun(i) = true;
                else
                    throw(ME);
                end
            end
        end
        
        function obj = remove_from_jobqueue(obj, ID, retry)
            % exact opposite of method add_from_jobqueue
            % Need to use JobID from obj.pool here instead of the jobqueue
            % index (obj.jobinfo.qi), which is not unique if uncomplete jobs exist from
            % previous modules
            
            ind = [obj.jobinfo.JobID] == ID;
            ji = obj.jobinfo(ind); % get job info struct
            
            % Backup Job Diary
            src = sprintf('%s/Job%d', obj.pool.JobStorageLocation, ji.JobID);
            dest = sprintf('%s_bck/Job%d', obj.pool.JobStorageLocation, ji.JobID);
            if exist(src,'dir')
                mkdir(dest);
                copyfile(src, dest);
            end
            
            % Clear job
            obj.jobinfo(ind) = [];
            obj.pool.Jobs([obj.pool.Jobs.ID] == ji.JobID).delete;
            
            % If retry requested, then reset jobnotrun and increment retry
            % counter
            if retry
                % Remove files from previous execution
                if exist(ji.jobpath,'dir'), rmdir(ji.jobpath,'s'); end
                
                % Add to jobqueue
                obj.jobretries.(ji.modulename)(ji.qi) = obj.jobretries.(ji.modulename)(ji.qi) + 1;
                obj.jobnotrun(ji.qi)=true;
            end
        end
                
        function states = job_monitor(obj, printjobs)
            % This function gathers job information from the job scheduler.
            % This can be slow depending on the size of the pool. 
            % INPUT
            % printjobs [true|false]: print job information to the screen.
            % OUTPUT
            % states: cell array of states (character) 
            %   {'running' | 'failed' | 'error' | 'finished'}
            % obj.jobinfo is also updated with the state information
            if nargin < 2
                printjobs = false;
            end
            if ~printjobs
                % This function might take a while. Let user know what's happening
                aas_log(obj.aap, false, 'Retrieving job information')
            end
            states = cell(1,numel(obj.jobinfo));
            jobids = [obj.pool.Jobs.ID];
            for id = shiftdim(jobids)'
                Jobs = obj.pool.Jobs(jobids == id);
                if isempty(Jobs) % cleared by the GUI
                    if obj.QV.isvalid
                        obj.QV.Hold = false;
                    end
                    obj.fatalerrors = true; % abnormal terminations
                    obj.close;
                    return;
                end
                jobind = [obj.jobinfo.JobID] == id;
                
                % If the job ID does not exist, jobinfo is not up-to-date
                % Assigning failed will cause the state handler to restart
                % the job without trying to remove it from pool.
                if any(jobind)
                    obj.jobinfo(jobind).state = Jobs.State;
                else
%                     obj.jobinfo(jobind).state = 'failed';
                end
                
                % Double check that finished jobs do not have an error in the Task object
                if strcmp(obj.jobinfo(jobind).state, 'finished')
                    if ~isempty(Jobs.Tasks.ErrorMessage)
                        obj.jobinfo(jobind).state = 'error';
                        
                    % Check if done flag exists. 
                    % Commented out. File system not always up to date and
                    % reliable.
%                     elseif ~exist(obj.jobqueue(obj.jobinfo(jobind).qi).doneflag, 'file')
%                         %                     obj.jobinfo(jobind).state = 'error';
%                         obj.jobinfo(jobind).state = 'error';
                    end
                end
            end
            
            states = {obj.jobinfo.state};

            if printjobs
                Nfinished = sum(strcmp(states,'finished'));
                Nqueued  = sum(strcmp(states,'queued'));
                Npending = sum(strcmp(states,'pending'));
                Nfailed = sum(strcmp(states,'failed'));
                Nerror    = sum(strcmp(states,'error'));
                Nrunning  = sum(strcmp(states,'running'));
                msg = sprintf('Running %3d | Queued %3d | Pending %3d | Finished %3d | Failed %3d | Error %3d',...
                    Nrunning, Nqueued, Npending, Nfinished, Nfailed, Nerror);
                aas_log(obj.aap,false,msg);
            end
            
            obj.QVUpdate;
        end
        
        function job_monitor_loop(obj)
            while true
                % Could put some backspaces here to remove the last
                % iteration
                states = obj.job_monitor(true); % states are also in e.g. obj.jobinfo(i).state
                if sum(strcmp(states, 'finished')) || sum(strcmp(states, 'error')) || sum(strcmp(states, 'failed'))
                    break;
                else
                    pause(10);
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UTILS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function esec = jts_etime(jts2,jts1)
    df = java.text.SimpleDateFormat('EEE MMM dd HH:mm:ss zzz yyyy');
    cal1 = java.util.GregorianCalendar.getInstance;
    cal2 = java.util.GregorianCalendar.getInstance;
    cal1.setTime(df.parse(jts1));
    cal2.setTime(df.parse(jts2));
    esec = etime([cal2.get(java.util.GregorianCalendar.YEAR) cal2.get(java.util.GregorianCalendar.MONTH)+1 cal2.get(java.util.GregorianCalendar.DAY_OF_MONTH) ...
        cal2.get(java.util.GregorianCalendar.HOUR_OF_DAY) cal2.get(java.util.GregorianCalendar.MINUTE) cal2.get(java.util.GregorianCalendar.SECOND)],...
        [cal1.get(java.util.GregorianCalendar.YEAR) cal1.get(java.util.GregorianCalendar.MONTH)+1 cal1.get(java.util.GregorianCalendar.DAY_OF_MONTH) ...
        cal1.get(java.util.GregorianCalendar.HOUR_OF_DAY) cal1.get(java.util.GregorianCalendar.MINUTE) cal1.get(java.util.GregorianCalendar.SECOND)]);
end