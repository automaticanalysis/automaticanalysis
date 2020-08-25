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
        function str = TaskInfo(obj,Task)
            modulename = Task.InputArguments{1}.tasklist.main.module(Task.InputArguments{3}).name;
            modality = Task.InputArguments{1}.schema.tasksettings.(modulename)(1).ATTRIBUTE.modality;
            acq = Task.InputArguments{1}.acq_details;
            indices = Task.InputArguments{4};
            
            switch modality
                case'MRI'
                    field_sess = 'sessions';
                    % for backwad compatibilty
                    if ~isempty(strfind(modulename,'diffusion'))
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
            str = sprintf('%s job %3d: %s',Task.State,Task.Parent.ID,Task.InputArguments{1}.tasklist.main.module(Task.InputArguments{3}).name);
        end
    end
end