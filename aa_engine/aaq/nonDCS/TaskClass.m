classdef TaskClass < handle
    properties
        Name
        CreateDateTime = datetime.empty
        StartDateTime = datetime.empty
        
        Parent
        
        InputArguments = {}        
    end
    
    properties (Hidden)
        Folder
        ShellFile
        LogFile
    end
    
    properties (Access = private)
        DiaryFile
        ProcessFile
        ErrorFile
    end
    
    properties (Dependent)
        State
        Worker
        FinishDateTime
        Diary
        ErrorMessage
        Error
    end
    
    methods
        function obj = TaskClass(Job,Name,varargin)
            obj.CreateDateTime = datetime('now','Timezone','local');
            
            obj.Parent = Job;
            obj.Name = sprintf('Task%d_%s',obj.Parent.latestTaskID+1,Name);
            obj.Folder = fullfile(obj.Parent.Folder,obj.Name);
            mkdir(obj.Folder)
            
            obj.ShellFile = fullfile(obj.Folder,'run.sh');
            obj.LogFile = fullfile(obj.Folder,'log.txt');
            obj.DiaryFile = fullfile(obj.Folder,'diary.txt');
            obj.ProcessFile = fullfile(obj.Folder,'process');
            obj.ErrorFile = fullfile(obj.Folder,'error.mat');
            
            % assemble command
            Command = func2str(varargin{1});
            if nargin < 4
                userVariable = [];                
            else
                obj.InputArguments = varargin{2};
                varStr = sprintf('arg%d,',1:numel(obj.InputArguments));
                varList = textscan(varStr,'%s','delimiter',',');
                userVariable = cell2struct(obj.InputArguments,varList{1}',2);                
                Command = [Command, '('];
                Command = [Command, varStr];
                Command(end) = ')';
            end                        
            if ~isempty(obj.Parent.AdditionalPaths)
                userVariable.reqpath = obj.Parent.AdditionalPaths;
                Command = sprintf('addpath(reqpath{:}); %s',Command);
            end            
            if ~isempty(userVariable)
                save(fullfile(obj.Folder,'data.mat'),'-struct','userVariable');
                Command = sprintf('load(''%s''); %s',fullfile(obj.Folder,'data.mat'),Command);
            end            
            Command(end+1) = ';';
            
            % create script
            fid = fopen(obj.ShellFile,'w');
            fprintf(fid,'export MALLOC_ARENA_MAX=4; matlab -nosplash -nodesktop -logfile %s -r "fid = fopen(''%s'',''w''); fprintf(fid,''%%s\\n'',java.lang.management.ManagementFactory.getRuntimeMXBean.getName.char); fclose(fid); try; %s catch E; save(''%s'',''E''); end; fid = fopen(''%s'',''a''); fprintf(fid,''%%s'',char(datetime(''now'',''Timezone'',''local''))); fclose(fid); quit"',...
                obj.DiaryFile,obj.ProcessFile,Command,obj.ErrorFile,obj.ProcessFile);
            fclose(fid);
        end
        
        function val = get.State(obj)
            val = 'unknown';
            if exist(obj.ProcessFile,'file')
                fid = fopen(obj.ProcessFile,'r');
                lines = textscan(fid,'%s','delimiter','@'); lines = lines{1};
                fclose(fid);
                if numel(lines) >= 3, val = 'finished'; end
            end
            if exist(obj.ErrorFile,'file'), val = 'error'; end
        end
        
        function val = get.Worker(obj)
            val = WorkerClass.empty;
            if exist(obj.ProcessFile,'file')
                fid = fopen(obj.ProcessFile,'r');
                lines = textscan(fid,'%s','delimiter','@'); lines = lines{1};
                fclose(fid);
                if numel(lines) >= 2, val = WorkerClass(lines{2},str2double(lines{1})); end
            end
        end
        
        function val = get.FinishDateTime(obj)
            val = datetime.empty;
            if any(strcmp({'finished','error'},obj.State))
                fid = fopen(obj.ProcessFile,'r');
                lines = textscan(fid,'%s','delimiter','@'); lines = lines{1};
                fclose(fid);
                if numel(lines) >= 3, val = datetime(lines{3},'Timezone','local'); end
            end
        end
        
        function val = get.Diary(obj)
            val = '';
            if exist(obj.DiaryFile,'file')
                fid = fopen(obj.DiaryFile,'r');
                while ~feof(fid)
                   val = cat(2,val,fgets(fid));
                end
                fclose(fid);
            end
        end
        
        function val = get.ErrorMessage(obj)
            val = '';
            if exist(obj.ErrorFile,'file')
                load(obj.ErrorFile)
                val = E.message;
            end
        end
        
        function val = get.Error(obj)
            val = MException.empty;
            if exist(obj.ErrorFile,'file')
                load(obj.ErrorFile)
                val = E;
            end
        end
        
    end
    
end

