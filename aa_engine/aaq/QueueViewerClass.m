% Queue Viewer Tool
% USAGE:
%   QV = QueueViewerClass(scheduler)
%       - scheduler: scheduler initiated (e.g. via parcluster)
%       - QV: queue viewer object
%
% GUI:
%   Clicking on jobs on the list provides info about the job
%   Clicking on "Info..." button provides info about the scheduler
%   Clicking on "Kill all!" button kills all jobs
%   Clicking on "Close" button closes the viewer. You can reopen using the Open method (e.g. QV.Open)
%
% Updating the list of jobs:
%   It does not automaticaly update the list of jobs after launching, but
%   you can update the list manually by clicking on the "Update" button or
%   calling the Update method (e.g. QV.Update). It also means that the the
%   program does not hold the focus, so you can continue working with MATLAB!
%
%   You can also switch on auto update by clicking on the "Auto update"
%   checkbox. With auto update swicthed on, you can set the update rate
%   using the slider. Keep in mind however, that in auto update mode, the
%   program holds the focus, so you cannot use MATLAB!
%
% Tibor Auer MRC CBU Cambridge 2012-2015

classdef QueueViewerClass < handle
    properties
        OnScreen
        Hold = true
    end
    properties (Access = protected)
        Pool
        DisplayHelper
        
        UIControls
                
        UpdateRate = 10
        Clock = []
    end
    properties (Access = private)
        JobsToIgnore = []
    end
    
    methods
        function obj = QueueViewerClass(pool)
            obj.Pool = pool;
            try
                obj.DisplayHelper = parallel.internal.display.DisplayHelper('Task');
            catch E
                if strcmp(E.identifier,'MATLAB:invalidType'), obj.DisplayHelper = parallel.internal.display.DisplayHelper(4);
                else rethrow(E); end
            end
            obj.Open;
        end
        
%         function delete(obj)
%             obj.Pool = [];
%         end
        
        function Open(obj)
            FontSize = [1.7*get(0,'DefaultUicontrolFontSize') get(0,'DefaultUicontrolFontSize')*0.8]; % H W
            DefaultListSize = [60 300]; % H W
            
            ncControls = 3; % #columns
            nrControls = 3; % #rows
            
            q = obj.GetQueue;
            
            [H, W] = size(char(q));
            
            % from listdlg
            fus = 8;
            ffs = 8;
            uh = 22; % Control Height
            
            figname = 'Queue';
            promptstring = {['Active jobs at ' datestr(now)]};
            listsize = [max([W*FontSize(2) DefaultListSize(2)]) max([H*FontSize(1) DefaultListSize(1)])];
            liststring = cellstr(q);
            initialvalue = 1;
            ex = FontSize(1);
            
            fp = get(0,'DefaultFigurePosition');
            w = 2*(fus+ffs)+listsize(1);
            h = 2*ffs+(5+nrControls)*fus+ex*numel(promptstring)+listsize(2)+nrControls*uh;
            fp = [fp(1) fp(2)+fp(4)-h w h];  % keep upper left corner fixed
            
            fig_props = { ...
                'name'                   figname ...
                'color'                  get(0,'DefaultUicontrolBackgroundColor') ...
                'resize'                 'off' ...
                'numbertitle'            'off' ...
                'menubar'                'none' ...
                'windowstyle'            'normal' ...
                'visible'                'off' ...
                'createfcn'              ''    ...
                'position'               fp   ...
                'closerequestfcn'        'delete(gcbf)' ...
                };
            obj.UIControls.fig = figure(fig_props{:});
            
            obj.UIControls.prompt_text = uicontrol('Style','text','String',promptstring,...
                'HorizontalAlignment','left',...
                'Position',[ffs+fus fp(4)-(ffs+fus+ex*numel(promptstring)) ...
                listsize(1) ex*numel(promptstring)]);
            
            obj.UIControls.listbox = uicontrol('Style','listbox',...
                'Position',[ffs+fus ffs+nrControls*uh+(3+nrControls)*fus listsize],...
                'String',liststring,...
                'BackgroundColor','w',...
                'Max',1,...
                'Tag','listbox',...
                'Value',initialvalue, ...
                'Callback', {@doListboxClick,obj});
            
            control_wid = (fp(3)-ncControls*(ffs+fus)-fus)/ncControls;
            
            pos = [1 1];
            obj.UIControls.autoupdate_lbl1 = uicontrol('Style','text','String','Autoupdate:',...
                'HorizontalAlignment','left',...
                'Position',[ffs+fus+(pos(2)-1)*(fus+control_wid) ffs+(nrControls-pos(1)+1)*fus+(nrControls-pos(1))*uh control_wid uh]);
            
            pos = [1 2];
            obj.UIControls.autoupdate_chk = uicontrol('Style','checkbox',...
                'Position',[ffs+fus+(pos(2)-1)*(fus+control_wid) ffs+(nrControls-pos(1)+1)*fus+(nrControls-pos(1))*uh control_wid uh],...
                'Tag','autoupdate_chk',...
                'Callback',{@chkAutoUpdate,obj});
            
            pos = [2 1];
            obj.UIControls.autoupdate_lbl2 = uicontrol('Style','text','String','Update rate:',...
                'Visible','off',...
                'HorizontalAlignment','left',...
                'Position',[ffs+fus+(pos(2)-1)*(fus+control_wid) ffs+(nrControls-pos(1)+1)*fus+(nrControls-pos(1))*uh control_wid uh]);
            
            pos = [2 2];
            obj.UIControls.autoupdate_sld = uicontrol('Style','slider',...
                'Visible','off',...
                'Min',10,'Max',24*3600,'Value',10,...
                'SliderStep',[10/(24*3600-10) 60/(24*3600-10)],...
                'Position',[ffs+fus+(pos(2)-1)*(fus+control_wid) ffs+(nrControls-pos(1)+1)*fus+(nrControls-pos(1))*uh control_wid*2 uh],...
                'Tag','autoupdate_sld',...
                'Callback',{@sldAutoUpdate,obj});
            
            pos = [1 3];
            obj.UIControls.update_btn = uicontrol('Style','pushbutton',...
                'String','Update',...
                'Position',[ffs+fus+(pos(2)-1)*(fus+control_wid) ffs+(nrControls-pos(1)+1)*fus+(nrControls-pos(1))*uh control_wid uh],...
                'Tag','update_btn',...
                'Callback',{@doUpdate,obj});
            
            pos = [3 1]; % row, column
            info_btn = uicontrol('Style','pushbutton',...
                'String','Info...',...
                'Position',[ffs+fus+(pos(2)-1)*(fus+control_wid) ffs+(nrControls-pos(1)+1)*fus+(nrControls-pos(1))*uh control_wid uh],...
                'Tag','info_btn',...
                'Callback',{@doInfo,obj});
            
            pos = [3 2];
            killall_btn = uicontrol('Style','pushbutton',...
                'String','Kill all!',...
                'Position',[ffs+fus+(pos(2)-1)*(fus+control_wid) ffs+(nrControls-pos(1)+1)*fus+(nrControls-pos(1))*uh control_wid uh],...
                'Tag','killall_btn',...
                'Callback',{@doKillAll,obj});
            
            pos = [3 3];
            close_btn = uicontrol('Style','pushbutton',...
                'String','Close',...
                'Position',[ffs+fus+(pos(2)-1)*(fus+control_wid) ffs+(nrControls-pos(1)+1)*fus+(nrControls-pos(1))*uh control_wid uh],...
                'Tag','close_btn',...
                'Callback',{@doClose,obj});
            
            set(obj.UIControls.fig, 'Visible','on'); drawnow;
            obj.OnScreen = true;
            obj.Hold = false;
        end
        
        function Close(obj)
            delete(obj.UIControls.fig);
            obj.OnScreen = false;
        end
        
        function Update(obj)
            q = obj.GetQueue;
            promptstring = get(obj.UIControls.prompt_text, 'String');
            promptstring{1} = ['Active jobs at ' datestr(now)];
            set(obj.UIControls.prompt_text, 'String',promptstring);
            set(obj.UIControls.listbox, 'String',cellstr(q));
            set(obj.UIControls.listbox, 'Value', min([1 numel(q)]))
            drawnow;
        end
        
        function UpdateAtRate(obj)
            if isempty(obj.Clock) || (toc(obj.Clock) >= obj.UpdateRate)
                obj.Clock = tic;
                obj.Update;
            end
        end
        
        function setAutoUpdate(obj,val)
            sw = {'off','on'};
            set(obj.UIControls.autoupdate_chk,'Enable',sw{val+1});
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UTILS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Hidden=true)
        function str = TaskInfo(obj,Task)
            elapsedtime = aas_getTaskDuration(Task);
            str = sprintf(['- Has been running on %s\n'...
                '- For %s'],...
                Task.Worker.Host,...
                elapsedtime);
        end
        
        function str = TaskLabel(obj,Task)
            str = sprintf('%s job %3d: %s',Task.State,Task.Parent.ID,func2str(Task.Function));
        end
        
        function q = GetQueue(obj)
            Jobs = [obj.Pool.Jobs];
            q = {};
            for j = 1:numel(Jobs)
                if any(obj.JobsToIgnore == Jobs(j).ID), continue; end
                Task = Jobs(j).Tasks;
                if ~strcmp(Task.State,'finished')
                    q{end+1} = obj.TaskLabel(Task);
                else % check for error
                    if ~isempty(Task.Error)
                        msg = sprintf('Job%d had an error: %s\n',Jobs(j).ID,Task.ErrorMessage);
                        for e = 1:numel(Task.Error.stack)
                            % Stop tracking to internal
                            if strfind(Task.Error.stack(e).file,'distcomp'), break, end
                            msg = [msg sprintf('in %s (line %d)\n', ...
                                Task.Error.stack(e).file, Task.Error.stack(e).line)];
                        end
                        msgbox(msg,[Jobs(j).Name Task.Name],'error')
                        if isfield(obj.UIControls, 'autoupdate_chk')
                            set(obj.UIControls.autoupdate_chk,'Value',false);
                            chkAutoUpdate(obj.UIControls.autoupdate_chk, [], obj)
                        end
                        obj.JobsToIgnore(end+1) = Jobs(j).ID;
                    end
                end
            end
            if isempty(q) && ~obj.Hold, msgbox('All jobs finished!','Queue','warn'); end
        end
        
        function KillAll(obj)
            for j = 1:numel(obj.Pool.Jobs)
                obj.Pool.Jobs(j).cancel;
            end
        end
        
        function str = TimeStr(obj,t)
            tdiv = [60 60 24 7];
            tformat = 'smhdw';
            str = '';
            
            for d = 1:numel(tdiv)
                td = rem(t,tdiv(d));
                if td, str = [num2str(td) tformat(d) ' ' str]; end
                t = floor(t/tdiv(d));
                if ~t, break; end
            end
        end
        
        function SetUpdateRate(obj,val)
            set(obj.UIControls.autoupdate_sld,'Value',val);
            obj.UpdateRate = get(obj.UIControls.autoupdate_sld,'Value');
            set(obj.UIControls.autoupdate_chk,'String',obj.TimeStr(obj.UpdateRate));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function chkAutoUpdate(autoupdate_chk, evd, obj) %#ok
switch get(autoupdate_chk,'Value')
    case 1 % auto
        set(obj.UIControls.update_btn,'Enable','off');
        set(obj.UIControls.autoupdate_lbl2,'Visible','on');
        set(obj.UIControls.autoupdate_sld,'Visible','on');
        set(obj.UIControls.fig, 'windowstyle','modal');
        set(autoupdate_chk,'String',obj.TimeStr(obj.UpdateRate));
        
        while get(autoupdate_chk,'Value')
            obj.UpdateAtRate;
            drawnow;
            if ~ishandle(obj.UIControls.fig), break; end
        end
    case 0 % manual
        set(obj.UIControls.update_btn,'Enable','on')
        set(obj.UIControls.autoupdate_lbl2,'Visible','off');
        set(obj.UIControls.autoupdate_sld,'Visible','off');
        set(obj.UIControls.fig, 'windowstyle','normal');
        set(autoupdate_chk,'String','');
end
end

function sldAutoUpdate(autoupdate_sld, evd, obj) %#ok
range = get(autoupdate_sld,'Max') - get(autoupdate_sld,'Min');

obj.SetUpdateRate(round(get(autoupdate_sld,'Value')/10)*10);
set(autoupdate_sld,'SliderStep',[10/range 60/range]);

if obj.UpdateRate > 60
    obj.SetUpdateRate(round(obj.UpdateRate/60)*60);
    set(autoupdate_sld,'SliderStep',[60/range 3600/range]);
end
end

function doUpdate(update_btn, evd, obj) %#ok
obj.Update;
end

function doKillAll(killall_btn, evd, obj) %#ok
obj.KillAll;
obj.Close;
delete(obj);
end

function doClose(close_btn, evd, obj) %#ok
obj.Close;
end

function doListboxClick(listbox, evd, obj) %#ok
% if this is a doubleclick, doOK
if strcmp(get(gcbf,'SelectionType'),'open')
    jobs = get(listbox,'String');
    sel = get(listbox,'Value');
    ID = regexp(jobs{sel},'job[ ]*[0-9]*:','match');
    ID = textscan(ID{1},'job %d:%*s'); ID = ID{1};
    Task = obj.Pool.Jobs([obj.Pool.Jobs.ID] == ID).Tasks;
    msgbox(obj.TaskInfo(Task),[Task.Parent.Name Task.Name],'help');
end
end

function doInfo(info_btn, evd, obj) %#ok
queue = 'unknown/unspecified';
mem = 'unknown/unspecified';
walltime = 'unknown/unspecified time';

switch class(obj.Pool)
    case 'parallel.cluster.Torque'
        if ~isempty(obj.Pool.SubmitArguments)
            schedinfo = textscan(obj.Pool.SubmitArguments,'%s'); schedinfo = schedinfo{1};
            if any(strcmp(schedinfo,'-q')), queue = schedinfo{circshift(strcmp(schedinfo,'-q'),1)}; end
            if cell_index(schedinfo,'mem='), mem = strrep(schedinfo{cell_index(schedinfo,'mem=')},'mem=',''); end
            if cell_index(schedinfo,'walltime='), walltime = obj.TimeStr(str2double(strrep(schedinfo{cell_index(schedinfo,'walltime=')},'walltime=',''))); end
        end
    case 'parallel.cluster.Generic'
        if any(strcmp(obj.Pool.IndependentSubmitFcn,'memory')), mem = sprintf('%d GB',obj.Pool.IndependentSubmitFcn{find(strcmp(obj.Pool.IndependentSubmitFcn,'memory'))+1}); end
        if any(strcmp(obj.Pool.IndependentSubmitFcn,'walltime')), walltime = sprintf('%d h',obj.Pool.IndependentSubmitFcn{find(strcmp(obj.Pool.IndependentSubmitFcn,'walltime'))+1}); end
end

msgbox(sprintf(['- Queuing %d jobs\n'...
    '- Running %d jobs\n'...
    '- Finished %d jobs\n'...
    '\n'...
    '- Submitted from %s\n'...
    '- To %d core(s)/Jobs\n'...
    '- In %s queue\n'...
    '- With %s RAM\n'...
    '- For maximum %s\n'...
    '\n'...
    '- Job Storage in %s'],...
    sum(strcmp({obj.Pool.Jobs.State},'queued')),...
    sum(strcmp({obj.Pool.Jobs.State},'running')),...
    sum(strcmp({obj.Pool.Jobs.State},'finished')),...
    obj.Pool.Host,...
    obj.Pool.NumWorkers,...
    queue,...
    mem,...
    walltime,...
    obj.Pool.JobStorageLocation),['Pool: ' class(obj.Pool)],'help');
end