classdef WorkerClass
    properties (SetAccess = immutable)
        Host
        ProcessId
    end
    
    methods
        function obj = WorkerClass(host,pid)
            obj.Host = host;
            obj.ProcessId = pid;
        end
    end
end

