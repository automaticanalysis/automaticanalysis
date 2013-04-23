function aap=aas_addinitialstream(aap,streamname,varargin)
%AAS_ADDINITIALSTREAM Add initial stream
%
% function aap=aas_addinitialstream(aap,streamname,varargin)
%
% If you use named subjects, you must put this command after the
% corresponding aas_addsubject command
%
% Example of use for named subject CBUtestme:
%
%  aap=aas_addinitialstream(aap,'structural','CBUtestme','/home/rcusack/my_structural.nii');
%
% Or, some other examples:
%
%   aap=aas_addinitialstream(aap,'structural',1,'/home/rcusack/Neonatal_fMRI/analaysis_v1/aamod_realign_00001/NNU1056/motor/meanNNU1056_fMRIbal_0001.nii');
%   aap=aas_addinitialstream(aap,'fieldmaps',3,{'/home/rcusack/my_fm1.nii','/home/rcusack/my_fm1.nii'});


% Set domain based on number of inputs
switch (length(varargin)-1)
    case 0
        domain='study';
    case 1
        domain='subject';
    case 2
        domain='session';
end;

if strcmp(domain, 'subject') || strcmp(domain, 'session')
    if ischar(varargin{1})
        subjnum=find(strcmp(varargin{1},{aap.acq_details.subjects.mriname}));
        if (isempty(subjnum))
            aas_log(aap,true,sprintf('Cannot find subject %s in list - is aas_addinitialstream command after aas_addsubject command?',varargin{1}));
        end;
    else
        subjnum=varargin{1};
    end;
end

if strcmp(domain, 'session')
    if ischar(varargin{2})
        sessnum=find(strcmp(varargin{2},{aap.acq_details.sessions.name}));
        if (isempty(sessnum))
            aas_log(aap,true,sprintf('Cannot find session %s in list - is aas_addinitialstream command after aas_addsession command?',varargin{1}));
        end;
    else
        sessnum=varargin{2};
    end;
end

% First, see if a module at the appropriate domain for this stream already
% exists

moduleexists=false;
for modposintasklist=1:length(aap.tasklist.main.module)
    if strcmp(aap.tasklist.main.module(modposintasklist).name,'aamod_importfilesasstream') ...
            && strcmp(streamname,aap.tasksettings.aamod_importfilesasstream(aap.tasklist.main.module(modposintasklist).index).outputstreams.stream) ...
            && strcmp(domain,aap.schema.tasksettings.aamod_importfilesasstream(aap.tasklist.main.module(modposintasklist).index).ATTRIBUTE.domain)
        moduleexists=true;
        break;
    end;
end;

if (~moduleexists)
    % CREATE NEW MODULE AND PUT AT START OF TASKLIST
    aap=aas_addtask(aap,'aamod_importfilesasstream',varargin(1:(end-1)));
    
    aap.aap_beforeuserchanges.tasksettings.aamod_importfilesasstream = aap.tasksettings.aamod_importfilesasstream;
    
    % Now we need to tailor with this
    if length(aap.tasklist.main.module)>1
        % First, lets move it to the beginning
        aap.tasklist.main.module=[aap.tasklist.main.module(end) aap.tasklist.main.module(1:(end-1))];
        % and do same to "before user changes" to prevent error trapping
        aap.aap_beforeuserchanges.tasklist.main.module=[aap.aap_beforeuserchanges.tasklist.main.module(end) aap.aap_beforeuserchanges.tasklist.main.module(1:(end-1))];
    end;
    
    if strcmp(domain, 'session')
        aap.tasklist.main.module(1).extraparameters.aap.acq_details.selected_sessions = sessnum;
        aap.aap_beforeuserchanges.tasklist.main.module(1).extraparameters.aap.acq_details.selected_sessions = sessnum;
    end
    
    % It's a little involved changing the engine's details by hand. For each
    % change, also make change to "aap_beforeuserchanges" to prevent error
    % trapping
    % Add the streams and the files to the task settings
    aap.tasksettings.aamod_importfilesasstream(end).outputstreams.stream{1}=streamname;
    aap.aap_beforeuserchanges.tasksettings.aamod_importfilesasstream(end).outputstreams.stream{1}=streamname;
    % And change the schema
    aap.schema.tasksettings.aamod_importfilesasstream(end).outputstreams.stream{1}=streamname;
    aap.aap_beforeuserchanges.schema.tasksettings.aamod_importfilesasstream(end).outputstreams.stream{1}=streamname;
    
    modposintasklist=1;
    modposinsettings=aap.tasklist.main.module(modposintasklist).index;
    aap.schema.tasksettings.aamod_importfilesasstream(modposinsettings).ATTRIBUTE.domain=domain;
else
    modposinsettings=aap.tasklist.main.module(modposintasklist).index;
end;


% ADD A MATCH ENTRY (SPECIFIC FOR THIS SUBJECT & SESSION IF NECESSARY)
% Add the files to import
fti=varargin{end};
if (~iscell(fti)), fti={fti}; end;
%matchind=length(aap.aap_beforeuserchanges.tasksettings.aamod_importfilesasstream(modposinsettings).match)+1;
matchind=length(aap.tasksettings.aamod_importfilesasstream(modposinsettings).match)+1;
aap.tasksettings.aamod_importfilesasstream(modposinsettings).match(matchind).filenames=fti;
aap.aap_beforeuserchanges.tasksettings.aamod_importfilesasstream(modposinsettings).match(matchind).filenames=fti;
% Set domain based on number of inputs
switch (length(varargin)-1)
    case 1
        aap.tasksettings.aamod_importfilesasstream(modposinsettings).match(matchind).subject=subjnum;
    case 2
        aap.tasksettings.aamod_importfilesasstream(modposinsettings).match(matchind).subject=subjnum;
        aap.tasksettings.aamod_importfilesasstream(modposinsettings).match(matchind).session=sessnum; %varargin{2};
end;
aap.aap_beforeuserchanges.tasksettings.aamod_importfilesasstream(modposinsettings).match=aap.tasksettings.aamod_importfilesasstream(modposinsettings).match;

