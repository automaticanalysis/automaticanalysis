classdef PoolClass < handle
    properties
        Type
        JobStorageLocation
        Jobs = JobClass.empty
        
        Host = getenv('HOSTNAME')
        Shell = 'bash'
        NumWorkers
        
        reqMemory = 1
        reqWalltime = 1
        initialConfiguration = ''
    end
    
    properties (Hidden)
        latestJobID = 0

        getSubmitStringFcn
        getSchedulerIDFcn
        getJobStateFcn
        getJobDeleteStringFcn
    end
    
    properties (Hidden, Access = protected)
        newGenericVersion
        initialSubmitArguments = ''
        SubmitArguments        
    end
    
    methods
        function obj = PoolClass(pool,initialSubmitArguments,initialConfiguration)
            % Required argument:
            % pool - Initializiation object (created by PCT function parcluster) or structure
            %
            % Required properties/fields of pool are:
            %
            % Type               - Scheduler 'Torque', 'LSF', 'Generic' (also for SGE)
            % JobStorageLocation - Location of log files for each job/task. Can be empty.
            %
            % Required properties/fields of pool for LSF:
            % SubmitArguments - Submission string pre-specifying memory and walltime. If empty defaults (1GB, 1hour) will be used.
            %
            % Required properties/fields of pool for Torque:
            % ResourceTemplate - Submission string pre-specifying memory and walltime. If empty defaults (1GB, 1hour) will be used.
            % SubmitArguments  - Submission string specifying resources other than memory and waltime. Can be empty.
            %
            % Required properties/fields of pool for Generic/SGE:
            % AdditionalProperties.AdditionalSubmitArgs - Submission string pre-specifying memory and walltime. If empty defaults (1GB, 1hour) will be used.
            %
            % Optionmal argument:
            % initialSubmitArguments - Submission string specifying resources other than memory and waltime. Can be empty.
            if nargin >= 2, obj.initialSubmitArguments = initialSubmitArguments; end
            if nargin >= 3, obj.initialConfiguration = initialConfiguration; end           
            if ~isempty(obj.initialSubmitArguments), obj.initialSubmitArguments = [' ' obj.initialSubmitArguments]; end
            
            obj.Type = pool.Type;
            obj.JobStorageLocation = pool.JobStorageLocation;
            obj.NumWorkers = pool.NumWorkers;
            
            switch obj.Type
                case 'Slurm'
                    obj.SubmitArguments = obj.initialSubmitArguments;
                    if isprop(pool,'ResourceTemplate')
                        obj.SubmitArguments = [obj.SubmitArguments ' ' pool.ResourceTemplate];
                    else
                        obj.SubmitArguments = [obj.SubmitArguments ' ' pool.SubmitArguments];
                    end
                    datWT = sscanf(regexp(obj.SubmitArguments,'-t [0-9]*','once','match'),'-t %d');
                    datMem = sscanf(regexp(obj.SubmitArguments,'--mem=[0-9]*[MGT]{1}','once','match'),'--mem=%d%c');
                    obj.getSubmitStringFcn = @(Job) sprintf( 'sbatch %s -J %s -o "%s" -e "%s" "%s"', ...
                         obj.SubmitArguments, Job.Name, Job.Tasks.LogFile, Job.Tasks.LogFile, Job.Tasks.ShellFile);
                    obj.getSchedulerIDFcn = @(stdOut) str2double(regexp(stdOut, '[0-9]*', 'match', 'once' ));
                    obj.getJobStateFcn = @(SchedulerID) Slurm_getJobState(SchedulerID);
                    obj.getJobDeleteStringFcn = @(SchedulerID) sprintf('scancel %d',SchedulerID);
                case 'Torque'
                    obj.SubmitArguments = obj.initialSubmitArguments;
                    if isprop(pool,'ResourceTemplate')
                        obj.SubmitArguments = pool.ResourceTemplate;
                    else
                        obj.SubmitArguments = pool.SubmitArguments;
                    end
                    datWT = sscanf(regexp(obj.SubmitArguments,'walltime=[0-9]*','once','match'),'walltime=%d');
                    datMem = sscanf(regexp(obj.SubmitArguments,'mem=[0-9]*','once','match'),'mem=%d');
                    obj.getSubmitStringFcn = @(Job) sprintf( 'qsub %s -N %s -j oe -o "%s" "%s"', ...
                        obj.SubmitArguments, Job.Name, Job.Tasks.LogFile, Job.Tasks.ShellFile);
                    obj.getSchedulerIDFcn = @(stdOut) str2double(regexp(stdOut, '[0-9]*', 'match', 'once' ));
                    obj.getJobStateFcn = @(SchedulerID) Torque_getJobState(SchedulerID);
                    obj.getJobDeleteStringFcn = @(SchedulerID) sprintf('qdel %d',SchedulerID);
                case 'LSF'
                    obj.SubmitArguments = pool.SubmitArguments;
                    datWT = sscanf(regexp(obj.SubmitArguments,'duration=[0-9]*','once','match'),'duration=%d');
                    datMem = sscanf(regexp(obj.SubmitArguments,'mem=[0-9]*','once','match'),'mem=%d')/1000;
                    obj.getSubmitStringFcn = @(Job) sprintf( 'bsub %s -J %s -oo "%s" < "%s"', ...
                        obj.SubmitArguments, Job.Name, Job.Tasks.LogFile, Job.Tasks.ShellFile);
                    obj.getSchedulerIDFcn = @(stdOut) str2double(regexp(stdOut, '[0-9]*', 'match', 'once' ));
                    obj.getJobStateFcn = @(SchedulerID) LSF_getJobState(SchedulerID);
                    obj.getJobDeleteStringFcn = @(SchedulerID) sprintf('bkill %d',SchedulerID);
                case 'Generic'
                    obj.newGenericVersion = ~isfield(pool,'IndependentSubmitFcn') || isempty(pool.IndependentSubmitFcn);
                    if obj.newGenericVersion
                        if ~isprop(pool.AdditionalProperties,'AdditionalSubmitArgs') && ~isfield(pool.AdditionalProperties,'AdditionalSubmitArgs')
                            warning(sprintf('WARNING: Propertiy "AdditionalSubmitArgs" not found.\n    "AdditionalSubmitArgs" must be listed within AdditionalProperties in the cluster profile in order to customise resource requirement and consequential queue selection.\n    Your jobs will be submitted to th default queue.'));
                        else
                            datWT = sscanf(regexp(pool.AdditionalProperties.AdditionalSubmitArgs,'s_cpu=[0-9]*','once','match'),'walltime=%d');
                            datMem = sscanf(regexp(pool.AdditionalProperties.AdditionalSubmitArgs,'s_rss=[0-9]*','once','match'),'mem=%d');
                        end
                    else
                        datWT = pool.IndependentSubmitFcn{find(strcmp(pool.IndependentSubmitFcn,'walltime'))+1};
                        datMem = pool.IndependentSubmitFcn{find(strcmp(pool.IndependentSubmitFcn,'memory'))+1};
                    end
                    obj.getSubmitStringFcn = @(Job) sprintf( 'qsub -S /bin/sh -N %s -j yes -o %s %s %s', ...
                        Job.Name, Job.Tasks.LogFile, obj.SubmitArguments, Job.Tasks.ShellFile);
                    obj.getSchedulerIDFcn = @(stdOut) sscanf(regexp(stdOut, 'Your job [0-9]*', 'once', 'match'),'Your job %d');
                    obj.getJobStateFcn = @(SchedulerID) SGE_getJobState(SchedulerID);
                    obj.getJobDeleteStringFcn = @(SchedulerID) sprintf('qdel %d',SchedulerID);
            end
            obj.reqWalltime = datWT;
            obj.reqMemory = datMem;
        end
        
        function set.JobStorageLocation(obj,value)
            obj.JobStorageLocation = value;
            if isempty(obj.JobStorageLocation)
                warning('JobStorageLocation is not specified. The current directory of %s will be used',pwd);
                obj.JobStorageLocation = pwd;
            elseif ~exist(obj.JobStorageLocation,'dir'), mkdir(obj.JobStorageLocation); 
            end
        end
            
        function set.reqWalltime(obj,value)
            if isempty(value), return; end
            obj.reqWalltime = value;
            obj.updateSubmitArguments;
        end
        
        function set.reqMemory(obj,value)
            if isempty(value), return; end
            obj.reqMemory = value;
            obj.updateSubmitArguments;
        end
        
        function Job = addJob(obj)
            Job = JobClass(obj);
            obj.Jobs(end+1) = Job;
            obj.latestJobID = obj.latestJobID + 1;
        end
        
        function val = get.Jobs(obj)
            val = obj.Jobs(obj.Jobs.isvalid);
        end
        
    end
    
    methods (Hidden, Access = protected)
        function updateSubmitArguments(obj)
            memory = obj.reqMemory;
            walltime = obj.reqWalltime;
            
            switch obj.Type
                case 'Slurm'
                    if round(memory) == memory % round
                        memory = sprintf('%dG',memory);
                    else % non-round --> MB
                        memory = sprintf('%dM',memory*1000);
                    end
                    obj.SubmitArguments = strcat(sprintf('--mem=%s -t %d ',memory,walltime*60),obj.initialSubmitArguments);
                case 'Torque'
                    if round(memory) == memory % round
                        memory = sprintf('%dGB',memory);
                    else % non-round --> MB
                        memory = sprintf('%dMB',memory*1000);
                    end
                    obj.SubmitArguments = strcat(sprintf('-q compute -l mem=%s -l walltime=%d',memory,walltime*3600),obj.initialSubmitArguments);
                case 'LSF'
                    obj.SubmitArguments = sprintf('%s -c %d -M %d -R "rusage[mem=%d:duration=%dh]"',obj.initialSubmitArguments,walltime*60,memory*1000,memory*1000,walltime);
                case 'Generic'
                    obj.SubmitArguments = sprintf('%s -l s_cpu=%d:00:00 -l s_rss=%dG',obj.initialSubmitArguments,walltime,memory);
            end
        end
    end
    
end

function state = SGE_getJobState(ID)
fname = tempname;
system(sprintf('qstat -xml > %s',fname));
qstat = xml_read(fname);

if isempty(qstat.queue_info), state = 'finished';return; end

joblist = qstat.queue_info.job_list;
ji = joblist([joblist.JB_job_number]==ID);
if isempty(ji)
    state = 'finished';
else
    state = ji.ATTRIBUTE.state; 
end
end

function state = Torque_getJobState(ID)
stateList = {...
    'HQ' 'queued';...
    'W' 'pending';...
    'R' 'running';...
    'C' 'finished';...
    'E' 'error'...
    };
[s, w] = system(sprintf('qstat -f %d',ID));
if s == 153 % Unknown Job ID
    state = 'finished';
else
    chState = regexp(w,'job_state = [A-Z]','match','once'); chState = chState(end); % RE HQW
    state = stateList{cell_index(stateList(:,1),chState),2};
end
end

function state = Slurm_getJobState(ID)
stateList = {...
    'PENDING' 'pending';...
    'RUNNING' 'running';...
    'COMPLETED' 'finished';...
    'CANCELLED' 'error';...
    'CANCELLED+' 'error';...
    'FAILED' 'error';...
    };
[s, w] = system(sprintf('sacct -j %d --format=state',ID)); w = textscan(w,'%s'); w = w{1};
if numel(w) == 2
    state = 'finished';
else
    state = stateList{cell_index(stateList(:,1),deblank(w{3})),2};
end
end

function state = LSF_getJobState(ID)
stateList = {...
    'PEND' 'pending';...
    'RUN' 'running';...
    'DONE' 'finished';...
    'EXIT' 'error'...
    };
[s, w] = system(sprintf('bjobs -noheader -o "stat" %d',ID));
if ~isempty(strfind(w,'not found'))
    state = 'finished';
else
    state = stateList{cell_index(stateList(:,1),deblank(w)),2};
end
end
