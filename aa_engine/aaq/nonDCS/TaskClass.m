classdef TaskClass < handle
    properties
        Name
        CreateTime
        FinishTime
        
        Parent
        
        InputArguments = {}        
        ShellFile
        LogFile
    end
    
    properties (Access = private)
        Folder
        DiaryFile
        DoneFile
        ErrorFile
    end
    
    properties (Dependent)
        State
        Diary
        ErrorMessage
        Error
    end
    
    methods
        function obj = TaskClass(Job,Name,varargin)
            obj.Parent = Job;
            obj.Name = sprintf('Task%d_%s',obj.Parent.latestTaskID+1,Name);
            obj.Folder = fullfile(obj.Parent.Folder,obj.Name);
            mkdir(obj.Folder)
            
            obj.ShellFile = fullfile(obj.Folder,'run.sh');
            obj.LogFile = fullfile(obj.Folder,'log.txt');
            obj.DiaryFile = fullfile(obj.Folder,'diary.txt');
            obj.DoneFile = fullfile(obj.Folder,'done');
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
            fprintf(fid,'matlab -nosplash -nodesktop -logfile %s -r "try; %s catch E; save(''%s'',''E''); end; fid = fopen(''%s'',''w''); fprintf(fid,''%%s'',char(toString(java.util.Date))); fclose(fid); quit"',...
                obj.DiaryFile,Command,obj.ErrorFile,obj.DoneFile);
            fclose(fid);
        end
        
        function val = get.State(obj)
            val = 'unknown';
            if exist(obj.DoneFile,'file')
                val = 'finished';
                fid = fopen(obj.DoneFile,'r');
                obj.FinishTime = fgetl(fid);
                fclose(fid);
            end
            if exist(obj.ErrorFile,'file'), val = 'error'; end
        end
        
        function val = get.Diary(obj)
            val = '';
            if exist(obj.DiaryFile,'file')
                fid = fopen(obj.DoneFile,'r');
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

