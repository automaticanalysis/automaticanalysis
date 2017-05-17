classdef aaq < handle
    properties
        aap
        isOpen      = false
        fatalerrors = false % flag to indicate fatal error --> stop pipeline
        jobqueue    = []
    end
    methods
        function obj=aaq(aap)
            if (exist('aap','var'))
                obj.aap=aap;
            end
            obj.isOpen = true;
        end
        
        function close(obj)
            aas_log(obj.aap,false,'Taskqueue is closed!');
            obj.isOpen = false;
        end
        
        %%==============================
        % Save self to a file
        function save(obj,fn)
            jobqueue=obj.jobqueue;
            aap=obj.aap;
            save(fn,'jobqueue','aap')
        end
        %%==============================
        % Clear the task queue
        function obj=emptyqueue(obj)
            obj.jobqueue=[];
        end
        
        %%==============================
        % Add a task to the task queue
        function obj=addtask(obj,taskmask)
            % Delete any to be completed firsts that have already been done
            tbcf={};
            for ind=1:length(taskmask.tobecompletedfirst)
                if (~aas_doneflagexists(obj.aap,taskmask.tobecompletedfirst{ind}))
                    tbcf{end+1}=taskmask.tobecompletedfirst{ind};
                end
            end
            taskmask.tobecompletedfirst=tbcf;
            
            obj.jobqueue=[obj.jobqueue,taskmask];
        end
        
        %%==============================
        % Allocate a job from the task queue to a worker
        function [obj couldbeallocated]=allocate(obj,i,highmem)
            global aaparallel;
            k=obj.jobqueue(i).k;
            [stagepath stagename]=fileparts(obj.aap.tasklist.main.module(k).name);
            try
                specialrequirements=obj.aap.tasksettings.(stagename)(obj.aap.tasklist.main.module(k).index).specialrequirements;
                %    specialrequirements={obj.aap.schema.tasksettings.(stagename)(obj.aap.tasklist.main.module(k).index).ATTRIBUTE.specialrequirements};
            catch
                specialrequirements={};
            end
            if exist('highmem','var') % djm: so can try with highmem after 1st crash, even if not set
                specialrequirements.highmemory=[];
                specialrequirements.unlimit=[];
            end
            
            
            [obj.aap, workerid]=aas_getboredworker(obj.aap,specialrequirements);
            
            couldbeallocated=false;
            if (~isempty(workerid))
                task=obj.getjobdescription(i);
                if (~isfield(aaparallel.workerstatus.(sprintf('worker%d',workerid)),'allocatedjobs'))
                    aaparallel.workerstatus.(sprintf('worker%d',workerid)).allocatedjobs={task};
                else
                    aaparallel.workerstatus.(sprintf('worker%d',workerid)).allocatedjobs=[aaparallel.workerstatus.(sprintf('worker%d',workerid)).allocatedjobs,{task}];
                end
                couldbeallocated=true;
            end
        end
        
        
        %% =================================
        % Get job description
        function [obj,task]=getjobdescription(obj,i)
            k=obj.jobqueue(i).k;
            clear task;
            task.name='doprocessing';
            task.aap=obj.aap;
            task.taskqueueposition=i;
            task.aap.tasklist.currenttask.epiprefix=obj.aap.tasklist.main.module(k).epiprefix;
            task.aap.tasklist.currenttask.extraparameters=obj.aap.tasklist.main.module(k).extraparameters;
            task.aap.internal.jobqueue=obj.jobqueue(i);
            task.aap.tasklist.currenttask.name=obj.aap.tasklist.main.module(k).name;
            task.aap.tasklist.currenttask.index=obj.aap.tasklist.main.module(k).index;
            % now set output root, which may be tailored for this stage
            task.aap.acq_details.root=obj.jobqueue(i).outputroot;
            task.aap.acq_details.rootsuffix=obj.jobqueue(i).rootsuffix;
        end
        
        %% The default, Mono threaded...
        
        % Run all tasks on the queue, single threaded
        function [obj]=runall(obj,dontcloseexistingworkers,waitforalljobs)
            global aaparallel
            
            njobs=length(obj.jobqueue);
            for i=1:njobs
                job=obj.jobqueue(i);
                obj.aap.acq_details.root=job.outputroot;
                obj.aap.acq_details.rootsuffix=job.rootsuffix;
                obj.aap.acq_details.inputrootsuffix=job.inputrootsuffix;
                
                aa_doprocessing_onetask(obj.aap,job.task,job.k);
            end
            obj.emptyqueue;
        end
       
    end
    
    %% Utils
    
    methods (Hidden, Access = protected)
        function argout = SetArg(obj,argin,key,value)
            argout = argin;
            if ~iscell(argout), argout = {argout}; end
            ind = find(cellfun(@(x) strcmp(x,key),argout));
            if ind, argout(ind:ind+1) = []; end
            argout(end+1:end+2) = {key value};
        end
    end
end
