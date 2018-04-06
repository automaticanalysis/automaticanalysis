classdef aaq_qsub_nonDCS < aaq_qsub
    
    methods
        function [obj]=aaq_qsub_nonDCS(aap)
            obj = obj@aaq_qsub(aap);
           
            if isempty(obj.pool), return; end
            
            pool = obj.pool;
            obj.pool = [];
            obj.pool = PoolClass(pool,obj.initialSubmitArguments);
            obj.pool.maximumRetry = obj.aap.options.aaworkermaximumretry;
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
                
                jobName = sprintf('aa%04d_%s',obj.pool.latestJobID+1,obj.aap.tasklist.main.module(job.k).name);
                % save job data
                aa.reqpath = aacache.path.reqpath;
                aa.aap = obj.aap;
                aa.job = job;
                aa.aaworker = aaworker;
                Job = obj.pool.addJob(jobName,...
                    'addpath(reqpath{:}); aa_doprocessing_onetask(aap,job.task,job.k,job.indices,aaworker);',...
                    aa);
                Job.Submit;
            else
                aa_doprocessing_onetask(obj.aap,job.task,job.k,job.indices);
            end
        end
        
    end
end