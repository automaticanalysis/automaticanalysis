classdef PoolClass < handle
    properties
        Type
        JobStorageLocation
        Jobs = JobClass.empty
        latestJobID = 0
        
        reqMemory = 1
        reqWalltime = 1
        
        getSubmitStringFcn
        getSchedulerIDFcn
        getJobStateFcn
        
        maximumRetry = 0
    end
    
    properties (Hidden, Access = protected)
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
                                        obj.SubmitArguments, Job.Name, Job.Log, Job.Script );
                    obj.getSchedulerIDFcn = @(stdOut) str2double(regexp(stdOut, '[0-9]*', 'match', 'once' ));
                    obj.getJobStateFcn = @(SchedulerID) PBS_getJobState(SchedulerID);
                case 'LSF'
                    obj.SubmitArguments = pool.SubmitArguments;
                case 'Generic'
                    obj.reqWalltime = pool.IndependentSubmitFcn{find(strcmp(pool.IndependentSubmitFcn,'walltime'))+1};
                    obj.reqMemory = pool.IndependentSubmitFcn{find(strcmp(pool.IndependentSubmitFcn,'memory'))+1};
                    obj.getSubmitStringFcn = @(Job) sprintf( 'qsub -S /bin/sh -N %s -j yes -o %s %s %s', ...
                        Job.Name, Job.Log, obj.SubmitArguments, Job.Script);
                    obj.getSchedulerIDFcn = @(stdOut) sscanf(regexp(stdOut, 'Your job [0-9]*', 'once', 'match'),'Your job %d');
                    obj.getJobStateFcn = @(SchedulerID) SGE_getJobState(SchedulerID);
            end
        end
        
        function set.reqWalltime(obj,value)
            obj.reqWalltime = value;
            obj.updateSubmitArguments;
        end
        
        function set.reqMemory(obj,value)
            obj.reqMemory = value;
            obj.updateSubmitArguments;
        end
        
        function Job = addJob(obj,jobName,Command,userVariable)
            if nargin < 4, userVariable = []; end
            Job = JobClass(obj,jobName,Command,userVariable);
            obj.Jobs(end+1) = Job;
            obj.latestJobID = obj.latestJobID + 1;
        end
        
        function submitJob(obj,Job)
            cmd = obj.SubmitString(Job);
            % submit
            % Job submission can sometimes fail (server fault) (DP). Added rety to cope this this.
            success = false;
            retries = 0;
            while success == false
                Job.Submit(cmd);
                success = ~strcmp(Job.State,'unknown');
                if ~success
                    if retries > obj.maximumRetry
                        error('#%d: %s',s,w);
                    else
                        warning('Could not add job. Retrying...')
                        retries = retries + 1;
                        pause(5)
                    end
                end
            end
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
chState = regexp(w,'job_state = [A-Z]','match','once'); chState = chState(end); % RE HQW
state = stateList{cell_index(stateList(:,1),chState),2};
end