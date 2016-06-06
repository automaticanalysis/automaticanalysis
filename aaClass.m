% aa version 5.2.0 May 2016
% Cusack R, Vicente-Grabovetsky A, Mitchell DJ, Wild CJ, Auer T, Linke AC, Peelle JE (2015). Automatic analysis (aa): Efficient neuroimaging workflows and parallel processing using Matlab and XML. Frontiers in Neuroinformatics 8:90.
% http://dx.doi.org/10.3389/fninf.2014.00090

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
        aaURL = 'http://automaticanalysis.org'
        aawiki = 'https://github.com/rhodricusack/automaticanalysis/wiki';
    end
    
    methods
        function obj = aaClass(varargin)
            obj.Name = 'aa';
            aafile = [mfilename('fullpath') '.m'];
            fid = fopen(aafile,'rt');
            if fid == -1, error('Can''t open %s.',aafile); end
            l1 = fgetl(fid); l2 = fgetl(fid); l3 = fgetl(fid);
            fclose(fid);
            d = textscan(l1,'%% aa version %s %s %d');
            obj.Version = d{1}{1};
            obj.Date = sprintf('%s %d',d{2}{1},d{3});
            obj.ManuscriptRef = l2(3:end);
            obj.ManuscriptURL = l3(3:end);
            obj.Path = fileparts(aafile);
            
            % Path
            if ~any(strcmp(varargin,'nopath'))
                fprintf('\nPlease wait a moment, adding <a href = "matlab: cd %s">%s</a> to the path\n',obj.Path,obj.Name);
                addpath(genpath(obj.Path)); % recursively add AA subfolders
                rmpath(genpath(fullfile(obj.Path,'.git'))); % remove GitHub-related path
            end
            
            % Greet
            if ~any(strcmp(varargin,'nogreet'))
                d = textscan(obj.ManuscriptRef,'%s %s %s','delimiter','.','CollectOutput',true); d = d{1};
                fprintf('Welcome to aa version %s %s\n',obj.Version,obj.Date);
                fprintf(' If you publish work that has used aa, please cite our manuscript:\n');
                fprintf(' <a href = "%s">%s</a>\n',obj.ManuscriptURL,d{1});
                fprintf(' <a href = "%s">%s</a>\n',obj.ManuscriptURL,d{2});
                fprintf(' <a href = "%s">%s</a>\n',obj.ManuscriptURL,d{3});
                fprintf('\nPlease visit <a href = "%s">The aa website</a> for more information!\n',obj.aaURL);
                fprintf('\nHere you can find example <a href = "matlab: cd %s">tasklists</a> and <a href = "matlab: cd %s">scripts</a>.\n',...
                    fullfile(obj.Path,'aa_recipes_and_parametersets'),fullfile(obj.Path,'examples'));
            end
            
            fprintf('Ready.\n');
        end
    end
end