classdef PoolClass < handle
    properties
        Type
        JobStorageLocation
        Jobs = JobClass.empty
        
        reqMemory = 1
        reqWalltime = 1
        
        maximumRetry = 0
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
        function obj = PoolClass(pool,initialSubmitArguments)
            if nargin >= 2, obj.initialSubmitArguments = initialSubmitArguments; end
            
            obj.Type = pool.Type;
            obj.JobStorageLocation = pool.JobStorageLocation;
            
            switch obj.Type
                case 'Torque'
                    obj.SubmitArguments = strcat(pool.ResourceTemplate, ' ', pool.SubmitArguments);
                    datWT = sscanf(regexp(obj.SubmitArguments,'walltime=[0-9]*','once','match'),'walltime=%d');
                    datMem = sscanf(regexp(obj.SubmitArguments,'mem=[0-9]*','once','match'),'mem=%d');
                    obj.reqWalltime = datWT;
                    obj.reqMemory = datMem;
                    obj.getSubmitStringFcn = @(Job) sprintf( 'qsub %s -N %s -j oe -o "%s" "%s"', ...
                                        obj.SubmitArguments, Job.Name, Job.Tasks.LogFile, Job.Tasks.ShellFile);
                    obj.getSchedulerIDFcn = @(stdOut) str2double(regexp(stdOut, '[0-9]*', 'match', 'once' ));
                    obj.getJobStateFcn = @(SchedulerID) PBS_getJobState(SchedulerID);
                    obj.getJobDeleteStringFcn = @(SchedulerID) sprintf('qdel %d',SchedulerID);
                case 'LSF'
                    obj.SubmitArguments = pool.SubmitArguments;
                case 'Generic'
                    obj.newGenericVersion = isempty(pool.IndependentSubmitFcn);
                    if obj.newGenericVersion
                        if ~isprop(pool.AdditionalProperties,'AdditionalSubmitArgs')
                            warning(sprintf('WARNING: Propertiy "AdditionalSubmitArgs" not found.\n    "AdditionalSubmitArgs" must be listed within AdditionalProperties in the cluster profile in order to customise resource requirement and consequential queue selection.\n    Your jobs will be submitted to th default queue.'));
                        else
                            datWT = sscanf(regexp(pool.AdditionalProperties.AdditionalSubmitArgs,'h_cpu=[0-9]*','once','match'),'walltime=%d');
                            datMem = sscanf(regexp(pool.AdditionalProperties.AdditionalSubmitArgs,'h_rss=[0-9]*','once','match'),'mem=%d');
                            obj.reqWalltime = datWT;
                            obj.reqMemory = datMem;
                        end
                    else
                        obj.reqWalltime = pool.IndependentSubmitFcn{find(strcmp(pool.IndependentSubmitFcn,'walltime'))+1};
                        obj.reqMemory = pool.IndependentSubmitFcn{find(strcmp(pool.IndependentSubmitFcn,'memory'))+1};
                    end
                    obj.getSubmitStringFcn = @(Job) sprintf( 'qsub -S /bin/sh -N %s -j yes -o %s %s %s', ...
                        Job.Name, Job.Tasks.LogFile, obj.SubmitArguments, Job.Tasks.ShellFile);
                    obj.getSchedulerIDFcn = @(stdOut) sscanf(regexp(stdOut, 'Your job [0-9]*', 'once', 'match'),'Your job %d');
                    obj.getJobStateFcn = @(SchedulerID) SGE_getJobState(SchedulerID);
                    obj.getJobDeleteStringFcn = @(SchedulerID) sprintf('qdel %d',SchedulerID);
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
                case 'Torque'
                    if round(memory) == memory % round
                        memory = sprintf('%dGB',memory);
                    else % non-round --> MB
                        memory = sprintf('%dMB',memory*1000);
                    end
                    obj.SubmitArguments = strcat(sprintf('-q compute -l mem=%s -l walltime=%d',memory,walltime*3600),obj.initialSubmitArguments);
                case 'Generic'
                    obj.SubmitArguments = sprintf('%s -l h_cpu=%d:00:00 -l h_rss=%dG',obj.initialSubmitArguments,walltime,memory);
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

function state = PBS_getJobState(ID)
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