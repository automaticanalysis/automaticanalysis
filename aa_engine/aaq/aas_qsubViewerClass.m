classdef aas_qsubViewerClass < QueueViewerClass
    methods
        function obj = aas_qsubViewerClass(varargin)
            if nargin
                taskqueue = varargin{1};
            else
                global taskqueue
            end
            obj = obj@QueueViewerClass(taskqueue.pool);
        end
    end
    
    methods (Hidden=true)
        function inargs = getInputArgs(~, Task)
            % Retrieve original input arguments to function 'createTask'
            % from input argument Task.
            % Background: function createTask is at the heart of the
            % aaq_qsub queue processor, and we need to retrieve the input
            % args to the function call here. createTask is either called
            % directly in aaq_qsub/qsub_q_job with the five input args
            % {obj.aap, job.task, job.k, job.indices, aaworker} or - in a
            % potential upcoming implementation using Matlab's batch
            % construct - from within BatchHelper2 via batch with a more
            % nested input arg structure.
            % We're identifying the implementation by looking for original
            % input arg aap, a struct with field 'directory_conventions'
            if isstruct(Task.InputArguments{1}) && isfield(Task.InputArguments{1}, 'directory_conventions')
                % original aaq_qsub implementation
                inargs = Task.InputArguments;
            elseif isstruct(Task.InputArguments{3}{1}) && isfield(Task.InputArguments{3}{1}, 'directory_conventions')
                % implementation using Matlab's 'batch'
                inargs = Task.InputArguments{3};
            end
        end
            
        function str = TaskInfo(obj,Task)
            inargs = obj.getInputArgs(Task);
            modulename = inargs{1}.tasklist.main.module(inargs{3}).name;
            modality = inargs{1}.schema.tasksettings.(modulename)(1).ATTRIBUTE.modality;
            acq = inargs{1}.acq_details;
            indices = inargs{4};
            
            switch modality
                case'MRI'
                    field_sess = 'sessions';
                    % for backwad compatibilty
                    if contains(modulename,'diffusion')
                        field_sess = 'diffusion_sessions';
                    end      
                case'DWI'
                    field_sess = 'diffusion_sessions';
                case { 'MEEG' 'MEG' 'EEG' }
                    field_sess = 'meeg_sessions';
            end
            if numel(indices) == 0, indicesstr = 'study';
            else
                if numel(indices) >= 1, indicesstr = sprintf('\n  - Subject %s',acq.subjects(indices(1)).subjname); end
                if numel(indices) >= 2, indicesstr = sprintf('%s\n  - Session %s',indicesstr,acq.(field_sess)(indices(2)).name); end
            end
            
            str = sprintf(['- Module: %s\n'...
                '- Indices: %s\n'...
                '%s'],...
                modulename,...
                indicesstr,...
                TaskInfo@QueueViewerClass(obj,Task));
        end
        
        function str = TaskLabel(obj,Task)
            inargs = obj.getInputArgs(Task);
            str = sprintf('%s job %3d: %s',Task.State,Task.Parent.ID,inargs{1}.tasklist.main.module(inargs{3}).name);
        end
    end
end