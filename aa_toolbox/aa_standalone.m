% Automatic analysis - Entry point for standalone version. 
%
% FORMAT function aa_standalone(fname_config, fname_tasklist, fname_aa [, 'mridatadir', <dir>][, 'meegdatadir', <dir>][, 'connection', <dir>][, 'anadir', <dir>][, 'subj', <subject indices>])
%   Compulsory arguments:
%       - fname_config: aa parameters (XML file)
%       - fname_tasklist: aa tasklist (XML file)
%       - fname_aa: pipeline customisation (XML file for editing aap structure)
%           - Its structure must correspond to that of the aap apart from special cases (see later).
%           - Index in lists of structures (e.g. aamod_smooth(1)) can be indicated following to aa practice: e.g. aamod_smooth_00001. Default is 1
%           - Special cases
%               - firstlevel_contrasts
%               - input.isBIDS
%   Optional arguments:
%       - 'mridatadir', <dir>: Directory to find raw MRI data (<dir> can be a colon separated list)
%       - 'meegdatadir', <dir>: Directory to find raw MEEG data
%       - 'connection', <dir>: Directory of a processed pipeline to connect to. If <dir> is a colon separated list (with two elements), then
%            first element is the directory of the the pipeline to connect to
%            second elements is the name of the latest stage to connect to
%       - 'anadir', <dir>: Directory to put analysis to
%       - 'subj', <subject name>: subjname(s) of subject(s) to be included

function aa_standalone(fname_config, fname_tasklist, fname_aa, varargin)

if strcmp(fname_config,'version'), aaClass('nopath'); return; end

%% Load tasklist and customisation
aap=aarecipe(fname_config,fname_tasklist);

%% Optiona arguments
args = vargParser(varargin);
if isfield(args,'mridatadir'), aap.directory_conventions.rawdatadir = args.mridatadir; end % data
if isfield(args,'meegdatadir'), aap.directory_conventions.rawmeegdatadir = args.meegdatadir; end % data
if isfield(args,'anadir') % output
    aap.acq_details.root = args.anadir; 
    aap.directory_conventions.analysisid = 'automaticanalysis'; 
end

%% User customisation
xml_aa = xml_read(fname_aa,struct('ReadAttr',0));
aap = recursive_set(aap,xml_aa);
    
% check path to T1 template --> the rest should work, too
spmdir = aas_gettoolboxdir(aap,'spm');
if exist(fullfile(spmdir,aap.directory_conventions.T1template),'file')
    aas_log(aap,false,['INFO: T1 template located in ' fullfile(spmdir,aap.directory_conventions.T1template)]);
else
    aas_log(aap,true,['T1 template cannot be found in ' fullfile(spmdir,aap.directory_conventions.T1template)]);
end
if isdeployed, aap.directory_conventions.templatedir = fullfile(ctfroot,aap.directory_conventions.templatedir); end
    
% BIDS
if xml_aa.acq_details.input.isBIDS
    aap = aas_processBIDS(aap); 
end

% constrast
fc = cellfun(@(x) regexp(x,'^firstlevel_contrasts.*','match'), fieldnames(xml_aa)','UniformOutput',false);
for mod = horzcat(fc{:})
    for con = fieldnames(xml_aa.(mod{1}))'
        if strcmp(con{1},'COMMENT') || ~isstruct(xml_aa.(mod{1}).(con{1})), continue; end
        aap = aas_addcontrast(aap,['aamod_' mod{1}],...
            xml_aa.(mod{1}).(con{1}).subject,...
            xml_aa.(mod{1}).(con{1}).format,...
            xml_aa.(mod{1}).(con{1}).vector,...
            xml_aa.(mod{1}).(con{1}).name,...
            xml_aa.(mod{1}).(con{1}).type);
    end
end

%% Subject selection
if isfield(args,'subj') 
    subjind = any(cell2mat(cellfun(@(x) strcmp({aap.acq_details.subjects.subjname},x)',cellstr(args.subj),'UniformOutput',false)),2);
    if ~any(subjind), aas_log(aap,false,'No subject found!'); end
    aap.acq_details.subjects = aap.acq_details.subjects(subjind);
end 

%% Connecting
if isfield(args,'connection') 
	con = textscan(args.connection,'%s','delimiter',':'); con = con{1};
	if numel(con) == 1, con{2} = ''; end % no maxstage specified
    
	aap=aas_doprocessing_initialisationmodules(aap);
    aap.directory_conventions.allowremotecache = 0;
    remotePipes = struct('host',           '', ...
        'directory',   con{1}, ...
        'allowcache',  0, ...
        'maxstagetag', con{2}, ...
        'checkMD5',    1);
    aap=aas_connectAApipelines(aap,remotePipes);
end

%% DO ANALYSIS
try
    aa_doprocessing(aap);
    aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
catch E
    aa_close(aap);
    rethrow(E);
end
end

function aap = recursive_set(aap,xml)
for f = fieldnames(xml)'
    if any(strcmp(f{1},{'firstlevel_contrasts','isBIDS','COMMENT'})), continue; end % special entries
    index = str2double(regexp(f{1},'[0-9]{5}','match')); if isempty(index), index = 1; end
    if isstruct(xml.(f{1})), aap.(f{1})(index) = recursive_set(aap.(f{1})(index),xml.(f{1})); 
    else
        aap.(f{1}) = xml.(f{1});
    end
end
end