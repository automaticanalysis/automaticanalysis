classdef JobClass < handle
    properties
        ID
        Name
        Folder
        Script
        Log
    end
    
    properties (Dependent)
        State
    end
    
    properties (Hidden, Access = protected)
        Pool
        schedulerID = NaN
    end
    
    methods
        function obj = JobClass(Pool,Name,Command,userVariable)
            obj.Pool = Pool;
            
            obj.Folder = fullfile(obj.Pool.JobStorageLocation,Name);
            while exist(obj.Folder,'dir')
                obj.Folder = [obj.Folder '+'];
            end
            mkdir(obj.Folder);            
            
            obj.ID = obj.Pool.latestJobID+1;
            obj.Name = Name;
            obj.Script = fullfile(obj.Folder,'run.sh');
            obj.Log = fullfile(obj.Folder,'log.txt');
            
            if (nargin >= 5) && ~isempty(userVariable)
                save(fullfile(obj.Folder,'data.mat'),'-struct','userVariable');
                Command = sprintf('load(''%s''); %s',fullfile(obj.Folder,'data.mat'),Command);
            end
            
            if Command(end) ~= ';', Command(end+1) = ';'; end
            
            % create script
            fid = fopen(obj.Script,'w');
            fprintf(fid,'matlab -nosplash -nodesktop -r "%s quit"',Command);
            fclose(fid);
        end
        
        function Submit(obj)
            [s, w] = system(obj.Pool.getSubmitStringFcn(obj));
            if ~s
                obj.schedulerID = obj.Pool.getSchedulerIDFcn(w);
            end            
        end
        
        function val = get.State(obj)
            val = 'unknown';
            if ~isnan(obj.schedulerID)
                val = obj.Pool.getJobStateFcn(obj.schedulerID);
            end           
        end
    end
    
end

