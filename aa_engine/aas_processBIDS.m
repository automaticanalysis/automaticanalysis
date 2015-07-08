% Allow specify Subjects, Sessions and Events from Brain Imaging Data Structure (BIDS)
% Required:
%     - aap.directory_conventions.rawdatadir: path to BIDS
% Optional:
%     - aap.acq_details.input.selected_sessions: selection of a subset of
%           tasks/runs based on their names (e.g. 'task001_run001')
%           (default = [], all tasks/runs are selected)
%     - aap.acq_details.input.correctEVfordummies: whether number of
%           dummies should be take into account when defining onset times (default = 1);
%           Also requires:
%               - aap.acq_details.numdummies: number of (acquired) dummy scans (default = 0);
%               - "repetition_time" (TR) specified in JSON header
%
% Tibor Auer MRC CBU Cambridge 2015

function aap = aas_processBIDS(aap)

%% Init
aap.directory_conventions.subjectoutputformat = '%s';
aap.directory_conventions.subject_directory_format = 3;
spm_jobman('initcfg');

if ~aap.directory_conventions.continueanalysis
    sfx = 1;
    analysisid0 = aap.directory_conventions.analysisid;
    while exist(fullfile(aap.acq_details.root,aap.directory_conventions.analysisid),'dir')
        sfx = sfx + 1;
        aap.directory_conventions.analysisid = [analysisid0 '_' num2str(sfx)];
    end
end

aap.tasksettings.aamod_firstlevel_model.xBF.UNITS  ='secs';

%% Process
% Set BIDS
BIDS = aap.directory_conventions.rawdatadir;

% Look for subjects
SUBJ = tsvread(fullfile(BIDS,'participants.tsv'));
iID = cell_index(SUBJ(1,:),'subject_id');

% Add
for subj = 2:size(SUBJ,1)
    subjID = SUBJ{subj,iID};
    SESS = spm_select('List',fullfile(BIDS,subjID),'dir','sess');
    if isempty(SESS)
        aap = add_data(aap,subjID,fullfile(BIDS,subjID));
    else
        for sess = 1:size(SESS,1)
            aap = add_data(aap,[subjID '_' SESS(sess,:)],fullfile(BIDS,subjID,SESS(sess,:)));
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UTILS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function aap = add_data(aap,subjstr,sesspath)

% Correct EV onstets for number of dummies?
numdummies = aap.acq_details.input.correctEVfordummies*aap.acq_details.numdummies;

images = {};
cf = cellstr(spm_select('List',sesspath,'dir'));
if any(strcmp(cf,'anatomy')), images = horzcat(images,cellstr(fullfile(sesspath,'anatomy',[subjstr '_T1w_001.nii.gz']))); end
if any(strcmp(cf,'functional'))
    for cfname = cellstr(spm_select('FPList',fullfile(sesspath,'functional'),[subjstr '_task.*_bold.nii.gz']))'
        taskname = strrep_multi(basename(cfname{1}),{[subjstr '_'] '_bold.nii'},{'',''});
        
        if ~isempty(aap.acq_details.input.selected_sessions) && ~any(strcmp(aap.acq_details.input.selected_sessions,taskname)), continue; end
        
        aap = aas_addsession(aap,taskname);
        
        % Data
        TR = 0;
        hdrfname = [strtok(taskname,'_') '_bold.json'];        
        hdrdir = sesspath;
        while ~exist(fullfile(hdrdir,hdrfname),'file') && ~strcmp(hdrdir,'/')
            hdrdir = fileparts(hdrdir);
        end
        if exist(fullfile(hdrdir,hdrfname),'file') 
            images = horzcat(images,struct('fname',cfname{1},'hdr',fullfile(hdrdir,hdrfname)));
            info = loadjson(fullfile(hdrdir,hdrfname));
            if isfield(info,'repetition_time'), TR = info.repetition_time; end
        else
            images = horzcat(images,cfname{1});
        end
        
        % Model
        if ~TR
            aas_log(aap,false,sprintf('WARNING: No header found for subject %s task/run %s\n',subjstr,taskname))
            aas_log(aap,false,'WARNING: No correction of EV onset for dummies is possible!')
        end

        % Search for event file
        eventfname = fullfile(sesspath,'functional',[subjstr '_' taskname '_events.tsv']); % default: next to the image
        if ~exist(eventfname,'file'), eventfname = fullfile(sesspath,'functional',[subjstr '_' strtok(taskname,'_') '_events.tsv']); end % same for all run
        if ~exist(eventfname,'file'), eventfname = fullfile(aap.directory_conventions.rawdatadir,[taskname '_events.tsv']); end % same for all subjects
        if ~exist(eventfname,'file'), eventfname = fullfile(aap.directory_conventions.rawdatadir,[strtok(taskname,'_') '_events.tsv']); end % same for all subjects and runs
        if ~exist(eventfname,'file')
            aas_log(aap,false,sprintf('WARNING: No event found for subject %s task/run %s\n',subjstr,taskname));
        else
            EVENTS = tsvread(eventfname);
            iName = cell_index(EVENTS(1,:),'trial_type');
            iOns = cell_index(EVENTS(1,:),'onset');
            iDur = cell_index(EVENTS(1,:),'duration');
            names = unique(EVENTS(2:end,iName)); onsets = cell(numel(names),1); durations = cell(numel(names),1);
            for t = 2:size(EVENTS,1)
                iEV = cell_index(names,EVENTS{t,iName});
                onsets{iEV}(end+1) = str2double(EVENTS{t,iOns});
                durations{iEV}(end+1) = str2double(EVENTS{t,iDur});
            end
            for e = 1:numel(names)
                aap = aas_addevent(aap,'aamod_firstlevel_model_00001',subjstr,taskname,names{e},onsets{e}-numdummies*TR,durations{e});
            end
        end;
    end
end
aap = aas_addsubject(aap,subjstr,images);
end

function LIST = tsvread(fname)
filestr = fileread(fname);
nLines = sum(double(filestr)==10);
LIST = textscan(filestr,'%s','Delimiter','\t'); 
LIST = reshape(LIST{1},[],nLines)';
end