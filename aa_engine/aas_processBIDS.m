% Allow specify Subjects, Sessions and Events from Brain Imaging Data Structure (BIDS)
% Required:
%     - aap.directory_conventions.rawdatadir: path to BIDS
% Optional:
%     - aap.acq_details.input.selected_sessions:
%           selection of a subset of tasks/runs based on their names (e.g. 'task001_run001')
%           - elements are strings: selecting for preprocessing and analysis
%           - elements are cells:
%               - one cell is a selection for one "aamod_firstlevel_model" (in the tasklist)
%               - tasks/runs selected for any "aamod_firstlevel_model" are pre-processed
%           (default = [], all tasks/runs are selected for both preprocessing and analysis)
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

for m = 1:numel(aap.tasksettings.aamod_firstlevel_model)
    aap.tasksettings.aamod_firstlevel_model(m).xBF.UNITS  ='secs';
end

%% Process
% Set BIDS
BIDS = aap.directory_conventions.rawdatadir;

% Look for subjects
SUBJ = tsvread(fullfile(BIDS,'participants.tsv'));
iID = cell_index(SUBJ(1,:),'subject_id');

% Add
for subj = 2:size(SUBJ,1)
    subjID = SUBJ{subj,iID};
    SESS = spm_select('List',fullfile(BIDS,subjID),'dir','ses');
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

% locate first_level modules
stagenumModel = cell_index({aap.tasklist.main.module.name},'aamod_firstlevel_model');

% Correct EV onstets for number of dummies?
numdummies = aap.acq_details.input.correctEVfordummies*aap.acq_details.numdummies;

images = {};
diffusionimages = {};
cf = cellstr(spm_select('List',sesspath,'dir'));
if any(strcmp(cf,'anatomy')), images = horzcat(images,cellstr(fullfile(sesspath,'anatomy',[subjstr '_T1w_001.nii.gz']))); end
if any(strcmp(cf,'diffusion'))
    for cfname = cellstr(spm_select('FPList',fullfile(sesspath,'diffusion'),[subjstr '_dwi_.*.nii.gz']))
        sessfname = strrep_multi(basename(cfname{1}),{[subjstr '_'] '.nii'},{'',''});
        sessname = sessfname;
        bvalfname = retrieve_file(fullfile(sesspath,'functional',[subjstr '_' sessfname '.bval']));
        bvecfname = retrieve_file(fullfile(sesspath,'functional',[subjstr '_' sessfname '.bvec']));
        if ~isempty(bvalfname) && ~isempty(bvecfname)
            aap = aas_add_diffusion_session(aap,sessname);
            diffusionimages = horzcat(diffusionimages,struct('fname',cfname{1},'bval',bvalfname,'bvec',bvecfname)); 
        else
            aas_log(aap,true,sprintf('ERROR: No BVals/BVecs found for subject %s run %s!\n',subjstr,sessfname))
        end
    end
end
if any(strcmp(cf,'functional'))
    for cfname = cellstr(spm_select('FPList',fullfile(sesspath,'functional'),[subjstr '_task.*_bold.nii.gz']))'
        taskfname = strrep_multi(basename(cfname{1}),{[subjstr '_'] '_bold.nii'},{'',''});
        info = []; TR = 0;
        taskname = taskfname;
        
        % Header
        hdrfname = retrieve_file(fullfile(sesspath,'functional',[subjstr '_' taskfname,'_bold.json']));
        if ~isempty(hdrfname)
            info = loadjson(hdrfname);
            if isfield(info,'TaskName')
                taskname = info.TaskName;
                taskname = strrep(taskname,' ',''); 
                taskname = strrep(taskname,',',''); 
                taskname = strrep(taskname,'(',''); 
                taskname = strrep(taskname,')',''); 
            end
            if isfield(info,'RepetitionTime'), TR = info.RepetitionTime; end
        end
        
        % Skip?
        if ~isempty(aap.acq_details.selected_sessions) && ~any(strcmp({aap.acq_details.sessions(aap.acq_details.selected_sessions).name},taskname)), continue; end
        aap = aas_addsession(aap,taskname);
        
        % Data
        if isstruct(info)
            images = horzcat(images,struct('fname',cfname{1},'hdr',hdrfname));
        else
            images = horzcat(images,cfname{1});
        end
        
        % Model
        if any(stagenumModel)
            iModel = [];
            for m = 1:numel(stagenumModel)
                sess = textscan(aap.tasklist.main.module(stagenumModel(m)).extraparameters.aap.acq_details.selected_sessions,'%s'); sess = sess{1};
                if any(strcmp(sess,'*')) || any(strcmp(sess,taskname)), iModel(end+1) = aap.tasklist.main.module(stagenumModel(m)).index; end
            end
            
            if ~TR
                aas_log(aap,false,sprintf('WARNING: No (RepetitionTime in) header found for subject %s task/run %s\n',subjstr,taskname))
                aas_log(aap,false,'WARNING: No correction of EV onset for dummies is possible!')
            end
            
            % Search for event file
            eventfname = retrieve_file(fullfile(sesspath,'functional',[subjstr '_' taskfname '_events.tsv'])); % default: next to the image
            if isempty(eventfname)
                aas_log(aap,false,sprintf('WARNING: No event found for subject %s task/run %s\n',subjstr,taskname));
            else
                EVENTS = tsvread(eventfname);
                iName = cell_index(EVENTS(1,:),'trial_type');
                iOns = cell_index(EVENTS(1,:),'onset');
                iDur = cell_index(EVENTS(1,:),'duration');
                names = unique(EVENTS(2:end,iName)); onsets = cell(numel(names),1); durations = cell(numel(names),1);
                for t = 2:size(EVENTS,1)
                    iEV = strcmp(names,EVENTS{t,iName});
                    onsets{iEV}(end+1) = str2double(EVENTS{t,iOns});
                    durations{iEV}(end+1) = str2double(EVENTS{t,iDur});
                end
                for m = iModel
                    for e = 1:numel(names)
                        aap = aas_addevent(aap,sprintf('aamod_firstlevel_model_%05d',m),subjstr,taskname,names{e},onsets{e}-numdummies*TR,durations{e});
                    end
                end
            end;
        end
    end
end
aap = aas_addsubject(aap,subjstr,images,[],[],diffusionimages);
end

function fname = retrieve_file(fname)
% fully specified fname with fullpath

% pre-process fname
[p, f, e] = fileparts(fname);
ntags = 1:(numel(find(fname=='_'))+1);
if isnan(str2double(f(end))) % last entry is not runnumber
    e = ['_' list_index(f,1,ntags(end)) e];
    ntags(end) = [];
    f = list_index(f,1,ntags);
end

% search
if ~exist(fname,'file'), fname = fullfile(fullfile(p,[list_index(f,1,ntags(1:end-1)) e])); end % try without runnumber
% one level up
p = fileparts(fileparts(p)); ntags(end-2) = [];
if ~exist(fname,'file'), fname = fullfile(fullfile(p,[list_index(f,1,ntags) e])); end
if ~exist(fname,'file'), fname = fullfile(fullfile(p,[list_index(f,1,ntags(1:end-1)) e])); end % try without runnumber
if numel(ntags) > 2 % we are in session folder
    % one level up
    p = fileparts(p); ntags(end-2) = [];
    if ~exist(fname,'file'), fname = fullfile(fullfile(p,[list_index(f,1,ntags) e])); end
    if ~exist(fname,'file'), fname = fullfile(fullfile(p,[list_index(f,1,ntags(1:end-1)) e])); end % try without runnumber
end
if ~exist(fname,'file'), fname = ''; end
end

function LIST = tsvread(fname)
filestr = fileread(fname);
nLines = sum(double(filestr)==10);
LIST = textscan(filestr,'%s','Delimiter','\t');
LIST = reshape(LIST{1},[],nLines)';
end

function str = strrep_multi(str, old, new)
for i = 1:numel(old)
    str = strrep(str, old{i}, new{i});
end
end