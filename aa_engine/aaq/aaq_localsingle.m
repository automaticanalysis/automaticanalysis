classdef aaq_localsingle<aaq        
    methods
        function [obj]=aaq_localsingle(aap)
            obj.aap=aap;
        end
%% The default, Mono threaded...
        
        % Run all tasks on the queue, single threaded
        function [obj]=runall(obj,dontcloseexistingworkers,waitforalljobs)
            global aaparallel
            
            njobs=length(obj.jobqueue);
            for i=1:njobs
                job=obj.jobqueue(i);
                aa_doprocessing_onetask(obj.aap,job.task,job.k,job.indices);
            end
            obj.emptyqueue;
        end
    end
end
        