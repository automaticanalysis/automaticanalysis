function varargout = userinput(varargin)
% Examples:
% resp = userinput('questdlg',sprintf('Cannot find parameters file %s\nSeed new parameter file from existing default?','paramfile'), 'Parameter file', 'Yes','No (Exit)','No (Exit)','GUI',true);
% [seedparam, rootpath] = userinput('uigetfile',{'*.xml','All Parameter Files' },'Desired seed parameter',defaultdir,'GUI',true);
% [defaultparameters, rootpath] = userinput('uiputfile',{'*.xml','All Parameter Files' }, 'Location of the parameters file and analyses by default',fullfile(pwd,defaultparameters),'GUI',true);

isGUI = true;
iParam = find(strcmpi(varargin,'gui'),1);
if ~isempty(iParam)
    isGUI = varargin{iParam+1};
    varargin(iParam:iParam+1) = [];
end

switch varargin{1}
    case 'questdlg' % question, title, btn1, btn2, (btn3,) default
        if isGUI
            mac_extra_print(varargin{3})
            varargout{1} = questdlg(varargin{2:end});
        else
            btns = varargin(4:end-1);
            msgBtn = sprintf(' / %s',btns{:}); msgBtn(1:3) = '';
            respList = cellfun(@(x) lower(strtok(x)), btns, 'UniformOutput', false);
            while true
                resp = input([varargin{2} ' (' msgBtn '):' ],'s');
                resp = btns(cellfun(@(x) strcmp(resp,x) || (resp==x(1)), respList));
                if ~isempty(resp), break; end
            end

            varargout{1} = resp{1};
        end
    case  'uigetdir' % path,title
        if isGUI
            mac_extra_print(varargin{3})
            varargout{1} = uigetdir(varargin{2:end});
        else
            % TODO: implement an abort option here too?
            rootpath = input([varargin{3} ' (or leave empty for ' varargin{2} '):'],'s');
            if isempty(rootpath), rootpath = varargin{2}; end

            varargout{1} = rootpath;
        end
    case  'uigetfile' % filter,title,defname
        if isGUI
            mac_extra_print(varargin{3})
            [varargout{1}, varargout{2}] = uigetfile(varargin{2:end});
        else
            defaultdir = varargin{4};
            defaultnames = dir(fullfile(defaultdir,varargin{2}{1}));
            fprintf('%s in %s:\n', varargin{2}{2}, defaultdir);
            fprintf('%s\n',defaultnames.name);
            while true
                seedparam = input([varargin{3} ' (or leave empty to abort):'],'s');
                % filter out extension (so we are robust to whether this is provided or not)
                [rootpath,seedparam,~] = fileparts(seedparam);
                if isempty(seedparam)
                    seedparam = 0; % 0 to indicate Cancel, same as uigetfile()
                    rootpath = 0;
                    break % leave empty to abort
                elseif isempty(rootpath)
                    rootpath = defaultdir;
                end
                seedparam = [seedparam varargin{2}{1}(2:end)]; %#ok<AGROW>

                if exist(fullfile(rootpath,seedparam),'file')
                    break
                else
                    fprintf('Could not find file %s in %s. Please try again!\n',seedparam, rootpath);
                end
            end

            varargout{1} = seedparam;
            varargout{2} = rootpath;
        end
    case  'uiputfile' % filter,title,defname
        if isGUI
            mac_extra_print(varargin{3})
            [varargout{1}, varargout{2}] = uiputfile(varargin{2:end});
        else
            % TODO: implement an abort option here too?
            [defaultdir, defaultseed]= fileparts(varargin{4});
            seedparam = input([varargin{3} ' (or leave empty for ' varargin{4} '):'],'s');
            % filter out extension (so we are robust to whether this is provided or not)
            [rootpath,seedparam] = fileparts(seedparam);
            if isempty(rootpath), rootpath = defaultdir; end
            if isempty(seedparam), seedparam = defaultseed; end

            varargout{1} = [seedparam varargin{2}{1}(2:end)];
            varargout{2} = rootpath;
        end
    otherwise
        error('Function %s is not an existing function or not implemented!',varargin{1});
end

end

function mac_extra_print(title)
% On mac, depending on os version, the custom dialog title does not always show
% Print the title to the command window before showing the dialog
if ismac()
    fprintf('\n');
    pause(0.5)
    fprintf('%s\n', title);
    pause(0.5)
end
end
