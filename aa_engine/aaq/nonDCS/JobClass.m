classdef JobClass < handle
    properties
        ID
        Name
        AdditionalPaths = {}
        Tasks = TaskClass.empty;
    end
    
    properties (Dependent)
        State
    end
    
    properties (Hidden)
        Folder
        latestTaskID = 0
    end
    
    properties (Hidden, Access = protected)
        Pool
        schedulerID = NaN
    end
    
    methods
        function obj = JobClass(Pool)
            obj.Pool = Pool;
            
            obj.ID = obj.Pool.latestJobID+1;
            obj.Name = sprintf('Job%d',obj.ID);
            
            obj.Folder = fullfile(obj.Pool.JobStorageLocation,obj.Name);
%             while exist(obj.Folder,'dir')
%                 obj.Folder = [obj.Folder '+'];
%             end
            mkdir(obj.Folder);            
        end
        
        function delete(obj)
            if ~obj.isvalid, return; end
            if obj.cancel
                rmdir(obj.Folder,'s'); 
                delete@handle(obj)
            else
                warning('Job %s could not be killed!\nYou may need to kill manually by calling %s.',obj.ID,obj.Pool.getJobDeleteStringFcn(obj.schedulerID));
            end
        end
        
        function s = cancel(obj)
            s = false;
            if ~any(strcmp({'unknown','finished','error'},obj.State)), [s, w] = system(obj.Pool.getJobDeleteStringFcn(obj.schedulerID)); end
            if s, warning(w); 
%             else
%                 obj.Pool.Jobs([obj.Pool.Jobs.ID]==obj.ID) = [];
            end
            s = ~s;
        end
        
        function addTask(obj,Name,varargin)
            Task = TaskClass(obj,Name,varargin{:});
            obj.Tasks(end+1) = Task;
            obj.latestTaskID = obj.latestTaskID + 1;
        end
        
        function Submit(obj)
            [s, w] = system(obj.Pool.getSubmitStringFcn(obj));
            obj.Tasks.StartDateTime = datetime('now','Timezone','local');
            if ~s
                obj.schedulerID = obj.Pool.getSchedulerIDFcn(w);
            end            
        end
        
        function val = get.State(obj)
            val = 'unknown';
            if ~isnan(obj.schedulerID)
                val = obj.Pool.getJobStateFcn(obj.schedulerID);
            end
            if any(strcmp({'finished','error'},val)), val = obj.Tasks.State; end
        end
    end
    
end

