classdef aaClass
    properties
        Path
    end
    
    properties (SetAccess = private)
        Name
        Version
        Date
        ManuscriptRef
        ManuscriptURL
        aaURL
        aawiki
    end
    
    methods
        function obj = aaClass(varargin)
            % Path
            aafile = [mfilename('fullpath') '.m'];
            obj.Path = fileparts(aafile);
            if ~any(strcmp(varargin,'nopath'))
                fprintf('\nPlease wait a moment, adding <a href = "matlab: cd %s">%s</a> to the path\n',obj.Path,obj.Name);
                addpath(genpath(obj.Path)); % recursively add AA subfolders
                rmpath(genpath(fullfile(obj.Path,'.git'))); % remove GitHub-related path
                rmpath(genpath(fullfile(obj.Path,'.github'))); % remove GitHub-related path
                rmpath(genpath(fullfile(obj.Path,'external','cluster'))); % remove cluster-integration path
                tbxdirs = dir(fullfile(obj.Path,'aa_tools','toolboxes'));
                % remove toolbox mods
                tbxdirs = tbxdirs(cellfun(@(d) ~isempty(regexp(d,'.*_mods$', 'once')), {tbxdirs.name}));
                rmpath(strrep(strjoin(cellfun(@genpath, fullfile(obj.Path,'aa_tools','toolboxes',{tbxdirs.name}), 'UniformOutput', false),pathsep),'::',':')); 
            end

            % user config directory
            configdir = fullfile(getenv('HOME'),'.aa');
            aas_makedir([], configdir);
            addpath(configdir);
            
            obj.Name = 'automaticanalysis';            
            info = loadjson(strrep(aafile,'Class.m','.json'));
            obj.Version = info.Version;
            obj.Date = info.Date;
            
            % get GitHub commit info if exists
            if exist(fullfile(obj.Path,'.git'),'dir')
                fid = fopen(fullfile(obj.Path,'.git','logs','HEAD'),'r');
                while ~(feof(fid)), line = fgetl(fid); end
                fclose(fid);
                % Split the line at tab (after tab = repo
                % location info)
                dat = textscan(line,'%s','delimiter','\t'); dat = dat{1};
                % Split at spaces (second to last var = date)
                dat = textscan(dat{1}, '%s','delimiter',' '); dat = dat{1};
                obj.Version = [obj.Version ' (' dat{2} ')'];
                obj.Date = datestr(str2double(dat{end-1})/86400 + datenum(1970,1,1),'mmm yyyy');
            end
            
            obj.ManuscriptRef = info.ManuscriptRef;
            obj.ManuscriptURL = info.ManuscriptURL;
            obj.aaURL = info.URL;
            obj.aawiki = info.wiki;            
            
            % Greet
            if ~any(strcmp(varargin,'nogreet'))
                d = textscan(obj.ManuscriptRef,'%s %s %s','delimiter','.','CollectOutput',true); d = d{1};
                fprintf('Welcome to aa version %s %s\n',obj.Version,obj.Date);
                fprintf(' If you publish work that has used aa, please cite our manuscript:\n');
                fprintf(' <a href = "%s">%s</a>\n',obj.ManuscriptURL,d{1});
                fprintf(' <a href = "%s">%s</a>\n',obj.ManuscriptURL,d{2});
                fprintf(' <a href = "%s">%s</a>\n',obj.ManuscriptURL,d{3});
                fprintf('\nPlease visit <a href = "%s">The aa website</a> for more information!\n',obj.aaURL);
                fprintf('\nHere you can find example <a href = "matlab: cd %s">parameter sets</a> and <a href = "matlab: cd %s">examples</a>.\n',...
                    fullfile(obj.Path,'aa_parametersets'),fullfile(obj.Path,'examples'));
                fprintf('Ready.\n');
            end
        end
    end
end
