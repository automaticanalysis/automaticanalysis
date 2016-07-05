% Allow specify Subjects, Sessions and Events from Brain Imaging Data Structure (BIDS)
% Required:
%     - aap.directory_conventions.rawdatadir: path to BIDS
% Optional:
%     - aap.acq_details.input.combinemultiple:
%           Combines multiple sessions per subjects (default = false)
%     - aap.acq_details.input.selected_sessions:
%           selection of a subset of tasks/runs based on their names (e.g. 'task-001_run01')
%           - elements are strings: selecting for preprocessing and analysis
%           - elements are cells:
%               - one cell is a selection for one "aamod_firstlevel_model" (in the tasklist)
%               - tasks/runs selected for any "aamod_firstlevel_model" are pre-processed
%           (default = [], all tasks/runs are selected for both preprocessing and analysis)
%     - aap.acq_details.input.correctEVfordummies: whether number of
%           dummies should be take into account when defining onset times (default = true);
%           Also requires:
%               - aap.acq_details.numdummies: number of (acquired) dummy scans (default = 0);
%               - "repetition_time" (TR) specified in JSON header
%
% Tibor Auer MRC CBU Cambridge 2015

function aap = aas_processBIDS(aap)

global BIDSsettings;
BIDSsettings.directories.structDIR = 'anat';
BIDSsettings.directories.functionalDIR = 'func';
BIDSsettings.directories.fieldmapDIR = 'fmap';
BIDSsettings.directories.diffusionDIR = 'dwi';
BIDSsettings.combinemultiple = aap.acq_details.input.combinemultiple;

%% Init
aap.directory_conventions.subjectoutputformat = '%s';
aap.directory_conventions.subject_directory_format = 3;

if ~aap.directory_conventions.continueanalysis
    sfx = 1;
    analysisid0 = aap.directory_conventions.analysisid;
    while exist(fullfile(aap.acq_details.root,aap.directory_conventions.analysisid),'dir')
        sfx = sfx + 1;
        aap.directory_conventions.analysisid = [analysisid0 '_' num2str(sfx)];
    end
end

if isfield(aap.tasksettings,'aamod_firstlevel_model')
    for m = 1:numel(aap.tasksettings.aamod_firstlevel_model)
        aap.tasksettings.aamod_firstlevel_model(m).xBF.UNITS  ='secs';
    end
end

%% Process
if BIDSsettings.combinemultiple
    aas_log(aap,false,'WARNING: You have selected combining multiple BIDS sessions!');
    aas_log(aap,false,'    Make sure that you have also set aap.options.autoidentify*_* appropriately!');
    aas_log(aap,false,'    N.B.: <aa sessionname> = <BIDS taskname>_<BIDS sessionname>');
end

% Set BIDS
BIDS = aap.directory_conventions.rawdatadir;

% Look for subjects
SUBJ = spm_select('List',aap.directory_conventions.rawdatadir,'dir','sub-.*');

% 1st pass - Add sessions only
% 2ns pass - Add data
for p = [false true]
    for subj = 1:size(SUBJ,1)
        subjID = deblank(SUBJ(subj,:));
        SESS = spm_select('List',fullfile(BIDS,subjID),'dir','ses');
        if isempty(SESS)
            aap = add_data(aap,subjID,fullfile(BIDS,subjID),p);
        else
            for sess = 1:size(SESS,1)
                aap = add_data(aap,[subjID '/' deblank(SESS(sess,:))],fullfile(BIDS,subjID,deblank(SESS(sess,:))),p);
            end
        end
    end
    if ~p && BIDSsettings.combinemultiple && (size(SESS,1) > 1)
        % sort sessions
        sessstr = regexp(cellstr(SESS),'-','split');
        sessstr = vertcat(sessstr{:});
        sessstr = sessstr(:,2);
        if numel(aap.acq_details.sessions) > 1
            aap.acq_details.sessions = aap.acq_details.sessions(sort_sessions(aap.acq_details.sessions,sessstr));
        end
        if numel(aap.acq_details.diffusion_sessions) > 1
            aap.acq_details.diffusion_sessions = aap.acq_details.diffusion_sessions(sort_sessions(aap.acq_details.diffusion_sessions,sessstr));
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UTILS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function aap = add_data(aap,mriname,sesspath,toAddData)

global BIDSsettings;
structDIR = BIDSsettings.directories.structDIR;
functionalDIR = BIDSsettings.directories.functionalDIR;
fieldmapDIR = BIDSsettings.directories.fieldmapDIR;
diffusionDIR = BIDSsettings.directories.diffusionDIR;

% locate first_level modules
stagenumModel = cell_index({aap.tasklist.main.module.name},'aamod_firstlevel_model');

% Correct EV onstets for number of dummies?
numdummies = aap.acq_details.input.correctEVfordummies*aap.acq_details.numdummies;

% Combine multiple sessions into one
if BIDSsettings.combinemultiple
    subjname = strtok(mriname,'/');
else
    subjname = strrep(mriname,'/','_');
end
aasubjname = strrep(subjname,'sub-','');

structuralimages = {};
functionalimages = {};
fieldmapimages = {};
diffusionimages = {};
for cf = cellstr(spm_select('List',sesspath,'dir'))'
    switch cf{1}
        case structDIR
            if ~toAddData, continue; end
            
            structuralimages = horzcat(structuralimages,cellstr(spm_select('FPList',fullfile(sesspath,structDIR),[subjname '.*_T1w.*.nii.gz'])));
        case diffusionDIR
            for cfname = cellstr(spm_select('FPList',fullfile(sesspath,diffusionDIR),[subjname '.*_dwi.*.nii.gz']))
                sessfname = strrep_multi(basename(cfname{1}),{[subjname '_'] '.nii'},{'',''});
                [bvalfname, runstr] = retrieve_file(fullfile(sesspath,diffusionDIR,[subjname '_' sessfname '.bval']));
                bvecfname = retrieve_file(fullfile(sesspath,diffusionDIR,[subjname '_' sessfname '.bvec']));
                if ~isempty(bvalfname) && ~isempty(bvecfname)
                    sessname = strrep_multi(sessfname,{basename(sesspath) runstr},{'',''}); sessname(1) = '';

                    if ~isempty(strfind(basename(sesspath),'ses-'))
                        sesstr = ['_' strrep(basename(sesspath),'ses-','')];
                    end
                    
                    if BIDSsettings.combinemultiple
                        sessname = [sessname,sesstr];
                    end
                    sessname = [sessname,runstr];
                    
                    aap = aas_add_diffusion_session(aap,sessname);
                    
                    if ~toAddData, continue; end
                    
                    aasessnames = {aap.acq_details.diffusion_sessions.name};
                    if BIDSsettings.combinemultiple
                        aasessnames = aasessnames(cell_index(aasessnames,sesstr));
                    end
                    if isempty(diffusionimages), diffusionimages = cell(1,numel(aasessnames)); end
                    isess = strcmp(aasessnames,sessname);
                    diffusionimages{isess} = struct('fname',cfname{1},'bval',bvalfname,'bvec',bvecfname);
                else
                    aas_log(aap,true,sprintf('ERROR: No BVals/BVecs found for subject %s run %s!\n',aasubjname,sessfname))
                end
            end
        case fieldmapDIR
            if ~toAddData, continue; end
            
            fmaps = cellstr(spm_select('FPList',fullfile(sesspath,fieldmapDIR),[subjname '.*.json']));
            skipnext = false;
            for f = fmaps'
                if skipnext
                    skipnext = false;
                    continue; 
                end
                jfname = spm_file(f{1},'basename');
                ind = find(jfname=='_',1,'last');
                ftype = jfname(ind+1:end);
                switch ftype
                    case 'phasediff'
                        fmap.hdr = loadjson(f{1});
                        if isfield(fmap.hdr,'IntendedFor')
                            fmap.hdr.session = get_taskname(sesspath,subjname,fmap.hdr.IntendedFor);
                        else
                            fmap.hdr.session = '*';
                        end
                        fmap.fname = cellstr(spm_select('FPList',fullfile(sesspath,fieldmapDIR),[strrep(jfname,ftype,'') '.*.nii.gz']));
                    case 'phase1'
                        skipnext = true;
                    case 'fieldmap'
                end
                fieldmapimages = horzcat(fieldmapimages,fmap);
            end            
        case functionalDIR
            for cfname = cellstr(spm_select('FPList',fullfile(sesspath,functionalDIR),[subjname '.*_task.*_bold.nii.gz']))'
                taskfname = strrep_multi(basename(cfname{1}),{[subjname '_'] '_bold.nii'},{'',''});
                [taskname, sesssfx] = get_taskname(sesspath,subjname,cfname{1});
                               
                % Header
                info = []; TR = 0;
                hdrfname = retrieve_file(fullfile(sesspath,functionalDIR,[subjname '_' taskfname,'_bold.json']));
                if ~isempty(hdrfname)
                    info = loadjson(hdrfname);
                    if isfield(info,'RepetitionTime'), TR = info.RepetitionTime; end
                end
                
                % Skip?
                if ~isempty(aap.acq_details.selected_sessions) && ~any(strcmp({aap.acq_details.sessions(aap.acq_details.selected_sessions).name},taskname)), continue; end
                aap = aas_addsession(aap,taskname);
                
                if ~toAddData, continue; end
                
                % Data
                aasessnames = {aap.acq_details.sessions.name};
                if BIDSsettings.combinemultiple
                    aasessnames = aasessnames(cell_index(aasessnames,sesssfx));
                end
                if isempty(functionalimages), functionalimages = cell(1,numel(aasessnames)); end
                isess = strcmp(aasessnames,taskname);
                if isstruct(info)
                    functionalimages{isess} = struct('fname',cfname{1},'hdr',hdrfname);
                else
                    functionalimages{isess} = cfname{1};
                end
                
                % Model
                if any(stagenumModel)
                    iModel = [];
                    for m = 1:numel(stagenumModel)
                        sess = textscan(aap.tasklist.main.module(stagenumModel(m)).extraparameters.aap.acq_details.selected_sessions,'%s'); sess = sess{1};
                        if any(strcmp(sess,'*')) || any(strcmp(sess,taskname)), iModel(end+1) = aap.tasklist.main.module(stagenumModel(m)).index; end
                    end
                    
                    if ~TR
                        aas_log(aap,false,sprintf('WARNING: No (RepetitionTime in) header found for subject %s task/run %s',aasubjname,taskname))
                        aas_log(aap,false,'WARNING: No correction of EV onset for dummies is possible!\n')
                    end
                    
                    % Search for event file
                    eventfname = retrieve_file(fullfile(sesspath,functionalDIR,[subjname '_' taskfname '_events.tsv'])); % default: next to the image
                    if isempty(eventfname)
                        aas_log(aap,false,sprintf('WARNING: No event found for subject %s task/run %s\n',aasubjname,taskname));
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
                                aap = aas_addevent(aap,sprintf('aamod_firstlevel_model_%05d',m),aasubjname,taskname,names{e},onsets{e}-numdummies*TR,durations{e});
                            end
                        end
                    end
                end
            end
        otherwise
            aas_log(aap,false,sprintf('NYI: Input %s is not supported',cf{1}));
    end
end
if toAddData
    aap = aas_addsubject(aap,aasubjname,mriname,...
        'structural',structuralimages,...
        'functional',functionalimages,...
        'fieldmaps',fieldmapimages,...
        'diffusion',diffusionimages); 
end
end

function [fname, runstr] = retrieve_file(fname)

% fully specified fname with fullpath

% pre-process fname
[p, f, e] = fileparts(fname);
tags = regexp(f,'_','split');
% suffix
if isempty(regexp(f(end),'[0-9]', 'once')) % last entry is suffix
    e = [tags{end} e];
    tags(end) = [];
end
ntags = 1:numel(tags);
f = list_index(f,1,ntags);
% suffix
itags = [cell_index(tags,'acq'), cell_index(tags,'run')]; itags(~itags) = [];
runstr = '';
for t = itags
    runstr = [runstr '_' tags{t}];
end

% search
itags = [cell_index(tags,'ses'), cell_index(tags,'sub')]; itags(~itags) = [];
itagrun = cell_index(tags,'run');
if ~exist(fname,'file') && itagrun, fname = fullfile(fullfile(p,[list_index(f,1,ntags(ntags~=itagrun),0) e])); end % try without runnumber
p = fileparts(p);
for t = itags
    % one level up
    p = fileparts(p); ntags(t) = [];
    if ~exist(fname,'file'), fname = fullfile(fullfile(p,[list_index(f,1,ntags,0) e])); end
    if ~exist(fname,'file') && itagrun, fname = fullfile(fullfile(p,[list_index(f,1,ntags(ntags~=itagrun),0) e])); end % try without runnumber
end
if ~exist(fname,'file'), fname = ''; end
end

function [taskname, sesstr] = get_taskname(sesspath,subjname,fname)
global BIDSsettings;
functionalDIR = BIDSsettings.directories.functionalDIR;
taskfname = strrep_multi(basename(fname),{[subjname '_'] '_bold.nii'},{'',''});
taskname = strrep(taskfname,'task-','');

% Header
[hdrfname, runstr] = retrieve_file(fullfile(sesspath,functionalDIR,[subjname '_' taskfname,'_bold.json']));
if ~isempty(hdrfname)
    info = loadjson(hdrfname);
    if isfield(info,'TaskName')
        taskname = info.TaskName;
        taskname = regexp(taskname,'[a-zA-Z0-9]*','match');
        taskname = strcat(taskname{:});
        if ~isempty(strfind(basename(sesspath),'ses-'))
            sesstr = ['_' strrep(basename(sesspath),'ses-','')];
        else
            sesstr = '';
        end
        if BIDSsettings.combinemultiple            
            taskname = [taskname,sesstr];
        end
        taskname = [taskname,runstr];
    end
end
end

function sessord = sort_sessions(sessions,sessstr)
aasessnames = {sessions.name};
sessstr = spm_file(sessstr,'prefix','_');
sessord = [];
for sess = sessstr
    sessord = horzcat(sessord, cell_index(aasessnames,sess{1}));
end
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