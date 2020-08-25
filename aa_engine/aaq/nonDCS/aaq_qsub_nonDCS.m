classdef aaq_qsub_nonDCS < aaq_qsub
    
    methods
        function [obj]=aaq_qsub_nonDCS(aap)
            global aaparallel;
            global aaworker;
            
            v0 = aap.options.verbose;
            aap.options.verbose = -1; % disable error
            obj = obj@aaq_qsub(aap);
            obj.aap.options.verbose = v0;

            if ~isempty(obj.pool)
                pool = obj.pool;
                obj.pool = [];
            else
                poolprofile = obj.poolConf{1};
                if exist(poolprofile,'file') % assume aap.directory_conventions.poolprofile is a path to settings
                    xml = xml_read(poolprofile);
                    SchedulerSettings = xml.settings(arrayfun(@(x) strcmp(x.ATTRIBUTE.name,'schedulercomponents'), xml.settings)).settings;
                    pool.Type = SchedulerSettings.ATTRIBUTE.name;
                    pool.JobStorageLocation =  aaworker.parmpath;
                    pool.ResourceTemplate = '';
                    pool.SubmitArguments = '';
                    pool.AdditionalProperties.AdditionalSubmitArgs = '';
                    pool.NumWorkers = SchedulerSettings.settings.key(arrayfun(@(x) strcmp(x.ATTRIBUTE.name,'NumWorkers'), SchedulerSettings.settings.key)).double.value;
                else
                    aas_log(obj.aap,false,'aap.directory_conventions.poolprofile is either not a file or cannot be found');
                    return
                end
            end
            obj.pool = PoolClass(pool,obj.initialSubmitArguments,obj.poolConf{3});
            obj.pool.Shell = aap.directory_conventions.linuxshell;
            obj.pool.reqMemory = aaparallel.memory;
            obj.pool.reqWalltime = aaparallel.walltime;
            obj.pool.NumWorkers = aaparallel.numberofworkers;
        end
        
        function obj = pool_args(obj,varargin)
            global aaparallel;
            memory = aaparallel.memory;
            walltime = aaparallel.walltime;
            
            for iarg = 1:numel(varargin)
                if ~ischar(varargin{iarg}), continue; end
                switch varargin{iarg}
                    case 'memory'
                        if ~isempty(varargin{iarg+1}), memory = varargin{iarg+1}; end
                    case 'walltime'
                        if ~isempty(varargin{iarg+1}), walltime = varargin{iarg+1}; end
                end
            end
            
            obj.pool.reqWalltime = walltime;
            obj.pool.reqMemory = memory;
        end
        
        
        function [obj]=qsub_q_job(obj,job)
            global aaworker
            global aacache
            aaworker.aacache = aacache;
            [s, reqpath] = aas_cache_get(obj.aap,'reqpath','system');
            % Let's store all our qsub thingies in one particular directory
            aas_makedir(obj.aap,fullfile(aaworker.parmpath,'qsub'));

            % Submit the job
            if ~isempty(obj.pool)
                % Check how much memory and time we should assign to the job
                qsubsettings = {'memory',[],'walltime',[]};
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
                obj = obj.pool_args(qsubsettings{:});
                
                taskName = obj.aap.tasklist.main.module(job.k).name;
                Job = obj.pool.addJob;
                Job.AdditionalPaths = reqpath;
                Job.addTask(taskName,@aa_doprocessing_onetask,{obj.aap,job.task,job.k,job.indices,aaworker});
                Job.Submit;
            else
                aa_doprocessing_onetask(obj.aap,job.task,job.k,job.indices);
            end
        end
        
    end
end