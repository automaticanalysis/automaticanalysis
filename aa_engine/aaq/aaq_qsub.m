classdef aaq_qsub<aaq
    properties
        pool = []
        QV = []
    end
    properties (Hidden)
        jobnotrun = []
        taskinqueue = []
        taskstomonitor = []
        jobinfo = []
        jobretries = []
        
        % ensure MAXFILTER license
        initialSubmitArguments = '';
    end
    methods
        function [obj]=aaq_qsub(aap)
            global aaworker;
            global aaparallel;
            aaparallel.numberofworkers=1;
            try
                if ~isempty(aap.directory_conventions.poolprofile)
                    profiles = parallel.clusterProfiles;
                    if ~any(strcmp(profiles,aap.directory_conventions.poolprofile))
                        ppfname = which(spm_file(aap.directory_conventions.poolprofile,'ext','.settings'));
                        if isempty(ppfname)
                            aas_log(obj.aap,true,sprintf('ERROR: settings for pool profile %s not found!',aap.directory_conventions.poolprofile));
                        else
                            obj.pool=parcluster(parallel.importProfile(ppfname));
                        end
                    else
                        aas_log(obj.aap,false,sprintf('INFO: pool profile %s found',aap.directory_conventions.poolprofile));
                        obj.pool=parcluster(aap.directory_conventions.poolprofile);
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
                        case 'parallel.cluster.Generic'
                            aas_log(obj.aap,false,'INFO: Generic engine is detected');
                            obj.pool.IndependentSubmitFcn = obj.SetArg(obj.pool.IndependentSubmitFcn,'walltime',aaparallel.walltime);
                            obj.pool.IndependentSubmitFcn = obj.SetArg(obj.pool.IndependentSubmitFcn,'memory',aaparallel.memory);
                    end
                else
                    obj.pool = parcluster('local');
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
            global aaworker
            
            % Check number of jobs & monitored files
            njobs=length(obj.jobqueue);
            
            % We have already submitted some of these jobs
            submittedJobs = 1:length(obj.jobnotrun);
            obj.jobnotrun = true(njobs,1);
            obj.jobnotrun(submittedJobs) = false;
            obj.jobretries = nan(njobs,1);
            jobqueuelimit = obj.aap.options.aaparallel.numberofworkers;
            printswitches.jobsinq = true; % switches for turning on and off messages
            
            
            while any(obj.jobnotrun) || (waitforalljobs && ~isempty(obj.jobinfo))
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
                            obj.job_monitor(false);
                        end
                    elseif ~isempty(obj.jobinfo)
                        % If no jobs left then monitor and update job states
                        aas_log(obj.aap, false, 'No jobs in the queue, waiting for remaining jobs to complete...')
                        obj.job_monitor_loop();
                        pause(60)
                    end
                else
                    aas_log(obj.aap, false, 'No free workers available: monitoring the queue...')
                    obj.job_monitor_loop()
                end
                
                idlist = [obj.jobinfo.ID];
                for id = idlist
                    ftmind = [obj.jobinfo.ID] == id;
                    state = obj.jobinfo(ftmind).state;
                    % This section of code only runs if the job is not
                    % running (has finished or error)
                    if strcmp(state, 'finished') || strcmp(state, 'error') || strcmp(state, 'failed')
                        Jobs = obj.pool.Jobs([obj.pool.Jobs.ID] == id);
                        if isempty(Jobs) % cleared by the GUI
                            if obj.QV.isvalid
                                obj.QV.Hold = false;
                            end
                            obj.fatalerrors = true; % abnormal terminations
                            obj.close;
                            return;
                        end
                        jobind = [obj.jobinfo.ID] == id;
                        %                     moduleName = obj.jobinfo(jobind).InputArguments{1}.tasklist.main.module(obj.jobinfo(jobind).InputArguments{3}).name;
                        aap = aas_setcurrenttask(obj.aap,obj.jobinfo(jobind).InputArguments{3});
                        moduleName = obj.aap.tasklist.main.module(obj.jobinfo(jobind).InputArguments{3}).name;
                        indices = obj.jobinfo(jobind).InputArguments{4};
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
                        
                        switch state
                            case 'pending'
                                t = toc(obj.jobinfo(jobind).tic);
                                % Tibor, you may want to add an option to
                                % aa to switch this on/off or extend the time? On very busy
                                % servers this might cause all jobs to be
                                % perpetually deleted and restarted.
                                if t > 3600 % if job has been pending for more than N seconds
                                    remove_from_jobqueue(obj, jobind)
                                end
                                
                            case 'failed' % failed to launch
                                msg = sprintf('Job%d had failed to launch (Licence?)!\n Check <a href="matlab: open(''%s'')">logfile</a>\n',id,...
                                    fullfile(obj.pool.JobStorageLocation,Jobs.Tasks.Parent.Name,[Jobs.Tasks.Name '.log']));
                                % If there is an error, it is fatal...
                                aas_log(obj.aap,false,msg,obj.aap.gui_controls.colours.error)
                                disp('Retrying')
                                remove_from_jobqueue(obj, jobind)
                                
                            case 'cancelled' % cancelled
                                msg = sprintf('Job%d had been cancelled by user!\n Check <a href="matlab: open(''%s'')">logfile</a>\n',id,...
                                    fullfile(obj.pool.JobStorageLocation,Jobs.Tasks.Parent.Name,[Jobs.Tasks.Name '.log']));
                                % If there is an error, it is fatal...
                                aas_log(obj.aap,true,msg,obj.aap.gui_controls.colours.warning)
                                
                            case 'finished' % without error
                                if isempty(Jobs.Tasks.FinishTime), continue; end
                                dtvs = dts2dtv(Jobs.Tasks.CreateTime);
                                dtvf = dts2dtv(Jobs.Tasks.FinishTime);
                                msg = sprintf('JOB %d: \tMODULE %s \tON %s \tSTARTED %s \tFINISHED %s \tUSED %s.',...
                                    id,moduleName,datname,Jobs.Tasks.CreateTime,Jobs.Tasks.FinishTime,sec2dts(etime(dtvf,dtvs)));
                                aas_log(obj.aap,false,msg,obj.aap.gui_controls.colours.completed);
                                
                                % Also save to file with module name attached!
                                fid = fopen(fullfile(aaworker.parmpath,'qsub','time_estimates.txt'), 'a');
                                fprintf(fid,'%s\n',msg);
                                fclose(fid);
                                
                                remove_from_jobqueue(obj, jobind, true)
                                
                            case 'error' % running error
                                
                                % Check whether the error was a "file does not
                                % exist" type. This can happen when a dependent
                                % folder is only partially written upon job execution.
                                % jobinfo etc is indexed by the job ID so get i from jobinfo.
                                
                                if obj.jobretries(jobind) <= 5
                                    obj.jobnotrun(jobind) = true; % will cause the job to restart on next loop
                                    msg = sprintf(['%s\n\n JOB FAILED WITH ERROR: \n %s',...
                                        ' \n\n Waiting 60 seconds then trying again',...
                                        ' (%d tries remaining for this job)\n'...
                                        'Press Ctrl+C now to quit, then run aaq_qsub_debug()',...
                                        ' to run the job locally in debug mode.\n'],...
                                        Jobs.Tasks.Diary, Jobs.Tasks.ErrorMessage, 5 - obj.jobretries(jobind));
                                    aas_log(aap, false, msg);
                                    pause(60)
                                    Jobs.delete; % delete this job from the cluster
                                    obj.jobinfo(jobind) = []; % remove inputarguments (otherwise this can get quite large)
                                else
                                    msg = sprintf('Job%d on <a href="matlab: cd(''%s'')">%s</a> had an error: %s\n',id,datpath,datname,Jobs.Tasks.ErrorMessage);
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
                                end
                        end
                    end
                end
                
                % Loop if we are still waiting for jobs to finish...
                if waitforalljobs
                    if ~any(obj.jobnotrun) && sum(strcmp({obj.jobinfo.state}, 'running')) == 0
                        waitforalljobs = false;
                    end
                end
                
                if obj.aap.options.aaworkerGUI
                    % queue viewer
                    % if ~isempty(obj.QV) && ~obj.QV.isvalid % killed
                    %     return
                    % end
                    if ~isempty(obj.pool)
                        if (isempty(obj.QV) || ~obj.QV.OnScreen) % closed
                            obj.QV = aas_qsubViewerClass(obj);
                            obj.QV.Hold = true;
                            obj.QV.setAutoUpdate(false);
                        else
                            obj.QV.UpdateAtRate;
                            if waitforalljobs, obj.QV.Hold = false; end
                        end
                    end
                end
            end
        end
        
        function obj = QVClose(obj)
            if ~isempty(obj.QV) && obj.QV.isvalid
                obj.QV.Close;
                obj.QV.delete;
                obj.QV = [];
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
                %                 if isfield(obj.aap.tasksettings.(job.stagename)(obj.aap.tasklist.main.module(job.k).index),'qsub')
                %                     qsub = obj.aap.tasksettings.(job.stagename)(obj.aap.tasklist.main.module(job.k).index).qsub;
                %                     for f = fieldnames(qsub)'
                %                         switch f{1}
                %                             case 'memoryBase'
                %                                 qsubsettings{2} = qsub.memoryBase;
                %                             case 'timeBase'
                %                                 qsubsettings{4} = qsub.timeBase;
                %                         end
                %                     end
                %                 end
                
                if isa(obj.pool,'parallel.cluster.Torque'), obj = obj.pool_args(qsubsettings{:}); end
                
                J = createJob(obj.pool);
                cj = @aa_doprocessing_onetask;
                nrtn = 0;
                inparg = {job.aap,job.task,job.k,job.indices, aaworker};
                
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
                retries = 0;
                
                % Job submission can sometimes fail (server fault) (DP). Added retry to cope this this.
                while retries < 5
                    try
                        J.submit;
                        success = true;
                    catch ME
                        aas_log(obj.aap, false, sprintf('WARNING: Could not add job due to following error:\n\n%s\n\nRetrying...', ME.message))
                        retries = retries + 1;
                        pause(5)
                        success = false;
                    end
                end
                
                if ~success
                    aas_log(obj.aap, true, sprintf('ERROR: Could not add job due to following error:\n\n%s', ME.message))
                end
            else
                aa_doprocessing_onetask(obj.aap,job.task,job.k,job.indices);
            end
        end
        
        
        function obj = add_from_jobqueue(obj, i)
            global aaworker
            % Add a job to the queue
            job=obj.jobqueue(i);
            job.aap.acq_details.root=aas_getstudypath(job.aap,job.k);
            % Run the job
            obj.qsub_q_job(job);
            % Create job info for referencing later
            % (also clean up done jobs to prevent IDs occuring twice)
            latestjobid = max([obj.pool.Jobs.ID]);
            if ~all(obj.jobnotrun(i)) % if any jobs have been run yet
                obj.jobinfo([obj.jobinfo.ID] == latestjobid) = []; % remove prev job with same ID
            end
            
            ji.InputArguments = {[],job.task,job.k,job.indices, aaworker};
            ji.ID = latestjobid;
            ji.i = i;
            ji.jobrunreported = false;
            ji.state = 'pending';
            ji.tic = tic;
            obj.jobinfo = [obj.jobinfo, ji];
            
            obj.jobnotrun(i) = false;
            moduleName = obj.aap.tasklist.main.module(ji.InputArguments{3}).name;
            aas_log(obj.aap, false, sprintf('Added job %s with ID: %d | Jobs submitted: %d',moduleName, latestjobid, length(obj.pool.Jobs)))
        end
        
        function obj = remove_from_jobqueue(obj, i, finished)
            % exact opposite of method add_from_jobqueue
            % Sometimes necessary if the job has failed to reach pending
            ID = obj.jobinfo(i).ID;
            obj.backup_job_diary(ID);
            obj.jobinfo(i) = [];
            obj.pool.Jobs([obj.pool.Jobs.ID] == ID).delete;
            if ~finished
                obj.jobnotrun(i)=true;
            end
        end
        
        function backup_job_diary(obj, jobid)
            src = sprintf('%s/Job%d', obj.pool.JobStorageLocation, jobid);
            dest = sprintf('%s_bck/Job%d', obj.pool.JobStorageLocation, jobid);
            if exist(src,'dir')
                mkdir(dest)
                copyfile(src, dest);
            end
        end
        
        function states = job_monitor(obj,printjobs)
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
                jobind = [obj.jobinfo.ID] == id;
                if length(find(jobind)) > 1
                    aas_log(obj.aap,true,'Job ID conflict: you should restart aa');
                end
                obj.jobinfo(jobind).state = Jobs.State;
            end
            states = {obj.jobinfo.state};

            if printjobs
                Nfinished = sum(strcmp(states,'finished'));
                Npending  = sum(strcmp(states,'pending')) + sum(strcmp(states,'queued'));
                Nerror    = sum(strcmp(states,'error')) + sum(strcmp(states,'failed'));
                Nrunning  = sum(strcmp(states,'running'));
                msg = sprintf('Running %3.1d | Pending %3.1d | Finished %3.1d | Error %3.1d', Nrunning, Npending, Nfinished, Nerror);
                aas_log(obj.aap,false,msg);
            end
        end
        
        function job_monitor_loop(obj)
            fw = false;
            while fw == false
                % Could put some backspaces here to remove the last
                % iteration
                states = obj.job_monitor(true); % states are also in e.g. obj.jobinfo(i).state
                if sum(strcmp(states, 'finished')) || sum(strcmp(states, 'error')) || sum(strcmp(states, 'failed'));
                    fw = true;
                else
                    pause(10)
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UTILS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dtv = dts2dtv(dts)
s = textscan(dts,'%s'); s = s{1}; s(5) = [];
s = strcat(s,{' ',' ',' ',' ',' '}'); s = [s{:}]; s = s(1:end-1);
dtformat = 'ddd mmm dd HH:MM:SS yyyy';
dtv = datevec(s,dtformat);
end

function dts = sec2dts(dt)
dt_str = {'s','m','h'};
dt_div = [60 60 24];

dts = '';
for i = 1:numel(dt_str)
    dts = [' ' num2str(mod(dt,dt_div(i))) dt_str{i} dts]; dt = floor(dt/dt_div(i));
    if ~dt, break, end
end
if dt
    dts = [num2str(dt) 'd' dts];
else
    dts = dts(2:end);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DLG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
