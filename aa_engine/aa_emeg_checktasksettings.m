function [aap subblock doit resp settings]=aa_emeg_checktasksettings(taskpath,args)
% Check task inputs and collect settings if required.
%
% If task is run from aa, args will be: {aap, job [,subject [,block]])
% and job is likely to be 'checkrequirements' or 'doit'
% Task can also be run directly, in which case args will be empty
% If task was run from aaMEG, args were {aap [,<subject|'setup'> [,block]])}
%
% If args{2}='doit', then fills any unspecified settings with defaults, then runs task
% If args{2}='checkrequirements' & ui=1 & this task hasn't yet been run,
% then prompts to confirm or respecify settings, then returns without running task.
% If args{2} is nonexistent then prompts to confirm or respecify settings,
% then runs task.
% Otherwise just returns.
%
% Danny Mitchell 09/07/08

%% prepare outputs (if aap is empty it will be made later)
subblock=[];
doit=0;
resp=''; % what does aa use this for?
settings=[];

%% get aap structure and prompt for userscript if not found
[taskpth task]=fileparts(taskpath);
try aap=args{1};
catch
    addpath /imaging/local/spm/spm5
    addpath /imaging/local/spm/spm5/cbu_updates/
    addpath /imaging/local/meg_misc/
    [aap ok]=spm_select(1,'any','Please select aap_parameters file or script to generate aap structure','',pwd,'^aa.*\.m.*$');
    if ~ok; return; end
    [aappth aapnam ext]=fileparts(aap);
    if strcmp(ext,'.m')
        addpath(aappth); run(aap); rmpath(aappth);
    end

    if ~isstruct(aap);
        try load(aap);
        catch error('\nFailed to load aap structure.\n')
        end
    end

    args{1}=aap; args{2}='check&doit'; aap.options.userinterface=1;
end

%% check job and decide whether to use ui
try ui=aap.options.userinterface; catch ui=1; end
if strcmp(args{2},'doit');
    ui=0;
elseif strcmp(args{2},'domain')
    aas_log(aap,1,'\n Called task to get domain, but might be easier to get this from aap.schema? \n')
elseif strcmp(args{2},'description');
    % what does aa use this for?
elseif strcmp(args{2},'checkrequirements')
    % checkrequirements runs for all domains before 'doit'
    % I'll use it once per task for checking/filling settings via ui
    % If a subject or block has already been run, we probably don't need/want to
    % change any settings.
    if ui
        try
            if args{3}>1; aap=args{1}; return; end
            if args{4}>1; aap=args{1}; return; end
        catch % no args>1 so proceed to confirm task settings
        end
    else
        aap=args{1};
        try subblock=args(3); catch; end
        try subblock=[subblock args(4)]; catch; end
        if ~strcmp(aap.schema.tasksettings.(task).ATTRIBUTE.domain,'internal') ...
                || ~isempty(subblock)
            return;
        end
    end
end

%% prepare interactive window
if ui
    set(0,'units','pixels'); % eeglab changes this to normalised
    F=spm_figure('findWin','Interactive');
    if isempty(findstr('aa',get(F,'Name')))
        delete(F)
        F=spm_figure('CreateWin','Interactive','aaMEG','on');
        set(F,'color',[(str2double(sprintf('%g','d'))-97)/26 ...
            (str2double(sprintf('%g','j'))-97)/26 ...
            (str2double(sprintf('%g','m'))-97)/26]);
    end
end

%% get default settings and task attributes from xml file
xmlfile=fullfile(taskpth,[task '.xml']);
if ~exist(xmlfile,'file');
    error('\nCould not find %s\n',[task '.xml']);
end

try xml=xml_read(xmlfile);
catch error('\nFound %s, but failed to read it.\n',[task '.xml']);
end

try
    defs=xml.tasksettings.(task);
    try defs=rmfield(defs,'ATTRIBUTE'); catch end
    settings=fieldnames(defs);
    taskattributes=xml.tasksettings.(task).ATTRIBUTE;
catch
    % fprintf('\n Debug %s\n',mfilename); keyboard
    error('\nLoaded %s, but could not determine settings.\n',[task '.xml']);
end

%% get appropriate index for repeated tasks, or set to 1
try index=aap.tasklist.currenttask.index;
catch index=1;
end

%% get subject (and block) if needed but not specified
% might need to change this...
if strcmp(taskattributes.domain,'subject') || strcmp(taskattributes.domain,'session')
    try subblock=args(3);
    catch
        if ui
            try 
                try subblock=spm_input('Please select subject:','+1','m',char(aap.acq_details.subjects.megname),1:length(aap.acq_details.subjects),1);
                catch subblock=spm_input('Please select subject:','+1','m',char(aap.acq_details.subjects.mriname),1:length(aap.acq_details.subjects),1);
                end
            catch % if aap not filled/available, select directory
                [aap.acq_details.subjects.megname ok]=spm_select(1,'dir','Please select subject directory');
                if ok; subblock={1}; end
            end
        end
        if isempty(subblock); error('\nFailed to determine subject for task %s\n',task); end
    end
elseif strcmp(taskattributes.domain,'internal') && length(args)>2
    subblock=args(3);
else subblock=[];
end
if strcmp(taskattributes.domain,'session')
    try subblock=[subblock args(4)];
    catch
        if ui
            try
                subblock={subblock,spm_input('Please select block:','+1','m',char(aap.acq_details.sessions.name),aap.acq_details.selected_sessions,1)};
            catch % if aap not filled/available, select directory
                [aap.acq_details.sessions.name ok]=spm_select(1,'dir','Please select session directory',aap.acq_details.subjects.megname);
                if ok; subblock(2)=1; end
            end
        end
        if length(subblock)~=2; error('\nFailed to determine session for task %s\n.',task); end
    end
end

%% optionally display current task settings and decide how to respecify
if ui
    try
        aa_emeg_displaytaskparams(aap.tasksettings.(task)(index),sprintf('==Current options for %s(%g):',task,index));
        source=spm_input(sprintf('Confirm settings for task: %s',task),'+1','m', ...
            'From user script (specify if empty)|From user script (default if empty)|Specify all|Use defaults', ...
            {'script+specify','script+defaults','specify','defaults'},1);
    catch
        source='specify';
    end
else source='script+defaults';
end

%% go through aap and deal with missing fields, or replace all with default
for f=1:length(settings) % according to xml

    if strcmp(settings{f},'specialrequirements'); continue; end

    % check for dependency and ignore as appropriate (assume settings can
    % only be dependent on choices earlier in the list)
    if isfield(defs.(settings{f}).ATTRIBUTE,'dependency')
        if ~isfield(defs.(settings{f}).ATTRIBUTE,'onlyif')
            error('\nValue of dependency not specified for %s.\n',settings{f})
        else
            if ~isfield(aap.tasksettings.(task)(index),defs.(settings{f}).ATTRIBUTE.dependency)
                % required setting not present, so skip this setting
                aap.tasksettings.(task)(index).(settings{f})='n/a'; continue
            end
            value=aap.tasksettings.(task)(index).(defs.(settings{f}).ATTRIBUTE.dependency);
            valid=defs.(settings{f}).ATTRIBUTE.onlyif;
            if ischar(valid);
                valid=textscan(valid,'%s','delimiter','|,;');
                valid=valid{1};
            end
            switch value
                %case eval(defs.(settings{f}).ATTRIBUTE.onlyif) % dependency met
                case valid
                otherwise % dependency not met
                    aap.tasksettings.(task)(index).(settings{f})='n/a';
                    continue
            end
        end
    end

    if ~isempty(findstr(char(source),'script')) && isfield(aap.tasksettings.(task)(index),settings{f})
        continue % because using aap and this field exists
    end

    if isempty(findstr(char(source),'specify'))
        % fill empty field with default
        aap.tasksettings.(task)(index).(settings{f})=defs.(settings{f}).CONTENT;
    else % prompt for value
        % determine input type required (defaults to 'e' for evaluated)
        type='e';
        switch defs.(settings{f}).ATTRIBUTE.ui
            case {'text'}; type='s';
            case {'int'}; type='i';
            case {'double'}; type='r';
            case {'file'}; type='select_file';
            case {'dir'}; type='select_dir';
            case {'yesno'}; type='b'; options={'no','yes'}; values=[0 1];
            otherwise
                options={}; values=[];
                if ~isempty(strfind(defs.(settings{f}).ATTRIBUTE.ui,'|int'))
                    type='bi1';
                    options=textscan(defs.(settings{f}).ATTRIBUTE.ui,'%s','delimiter','|');
                    options=options{1}(1:end-1);
                elseif ~isempty(strfind(defs.(settings{f}).ATTRIBUTE.ui,'|'))
                    type='b';
                    options=textscan(defs.(settings{f}).ATTRIBUTE.ui,'%s','delimiter','|');
                    options=options{1};
                elseif ~isempty(strfind(defs.(settings{f}).ATTRIBUTE.ui,'_'))
                    type='m';
                    options=textscan(defs.(settings{f}).ATTRIBUTE.ui,'%s','delimiter','_');
                    options=options{1};
                end
        end
        switch type
            case 'select_file'
                defdir=fileparts(defs.(settings{f}).CONTENT);
                aap.tasksettings.(task)(index).(settings{f})=spm_select(1,'any',defs.(settings{f}).ATTRIBUTE.desc,defs.(settings{f}).CONTENT,defdir);
            case 'select_dir'
                aap.tasksettings.(task)(index).(settings{f})=spm_select(1,'dir',defs.(settings{f}).ATTRIBUTE.desc,defs.(settings{f}).CONTENT,defdir);
            case {'b', 'm', 'bi1'}
                if ~strcmp(type,'bi1') && isempty(values)
                    try
                        values=textscan(defs.(settings{f}).ATTRIBUTE.values,'%s','delimiter','|,');
                        values=values{1};
                    catch values=[];
                    end
                    if length(values)~=length(options); values=1:length(options); end;
                end
                if isinf(defs.(settings{f}).CONTENT); default=Inf;
                else default=uint8(defs.(settings{f}).CONTENT);
                end
                if ~all(size(default)==1) | default<1; default=1; end % vector
                choice=spm_input(defs.(settings{f}).ATTRIBUTE.desc,'+1',type,options,values,default);
                if iscell(choice); choice=choice{1}; end
                aap.tasksettings.(task)(index).(settings{f})=choice;
            otherwise
                aap.tasksettings.(task)(index).(settings{f})=spm_input(defs.(settings{f}).ATTRIBUTE.desc,'+1',type,defs.(settings{f}).CONTENT);
        end
    end
end

%% update currenttask (?);
aap.tasklist.currenttask.settings=aap.tasksettings.(task)(index);
% aap.tasklist.currenttask.index=index; aap.tasklist.currenttask.name=task;

%% return settings and doit flag; report any action
settings=aap.tasklist.currenttask.settings;
if ~isempty(strfind(args{2},'doit'))
    if ~strcmp('internal',taskattributes.domain);
        aas_emeg_reportcurrentjob(aap,task,subblock);
    end
    doit=1; dbstop if error
elseif ~isempty(strfind(args{2},'parallelise'));
    aas_emeg_reportcurrentjob(aap,task,subblock,'Parallelising');
    doit=1; dbstop if error
else doit=0;
end

%% save aap ?

return