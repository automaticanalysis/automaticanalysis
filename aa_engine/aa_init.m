% Automatic analysis - initialise paths from recipe

function [aap]=aa_init(aap)

global aacache
aacache.bcp_path = path;

% Get root aa directory
mfp=fileparts(fileparts(which('aa_doprocessing')));

%% Set Paths
% Path for SPM
if isempty(aap.directory_conventions.spmdir)
    if isempty(which('spm'))
        aas_log(aap,true,'You''re going to need SPM, add it to your paths manually or set aap.directory_conventions.spmdir');
    else
        aap.directory_conventions.spmdir=spm('Dir');
    end;
end;

addpath(aap.directory_conventions.spmdir);
    
spm_jobman('initcfg');

if ~isfield(aap,'spm') || ~isfield(aap.spm,'defaults')
    try
        aap.spm.defaults=spm_get_defaults;
    catch
        global defaults
        if (~isstruct(defaults))
            aas_log(aap,true,'Global SPM defaults has not been found;');
        else
            aap.spm.defaults=defaults;
        end;
    end;
    
    % Make copy of aap
    aap.aap_beforeuserchanges.spm.defaults = aap.spm.defaults;
end

% Path for SPM MEG/EEG
addpath(fullfile(spm('Dir'),'external','fieldtrip'));
clear ft_defaults
clear global ft_default
ft_defaults;
global ft_default
ft_default.trackcallinfo = 'no';
ft_default.showcallinfo = 'no';
addpath(...
    fullfile(spm('Dir'),'external','bemcp'),...
    fullfile(spm('Dir'),'external','ctf'),...
    fullfile(spm('Dir'),'external','eeprobe'),...
    fullfile(spm('Dir'),'external','mne'),...
    fullfile(spm('Dir'),'external','yokogawa_meg_reader'),...
    fullfile(spm('Dir'),'toolbox', 'dcm_meeg'),...
    fullfile(spm('Dir'),'toolbox', 'spectral'),...
    fullfile(spm('Dir'),'toolbox', 'Neural_Models'),...
    fullfile(spm('Dir'),'toolbox', 'MEEGtools'));

% Path for EEGLAB, if specified
if ~isempty(aap.directory_conventions.eeglabdir)
    addpath(...
        fullfile(aap.directory_conventions.eeglabdir,'functions'),...
        fullfile(aap.directory_conventions.eeglabdir,'functions', 'adminfunc'),...
        fullfile(aap.directory_conventions.eeglabdir,'functions', 'sigprocfunc'),...
        fullfile(aap.directory_conventions.eeglabdir,'functions', 'guifunc'),...
        fullfile(aap.directory_conventions.eeglabdir,'functions', 'studyfunc'),...
        fullfile(aap.directory_conventions.eeglabdir,'functions', 'popfunc'),...
        fullfile(aap.directory_conventions.eeglabdir,'functions', 'statistics'),...
        fullfile(aap.directory_conventions.eeglabdir,'functions', 'timefreqfunc'),...
        fullfile(aap.directory_conventions.eeglabdir,'functions', 'miscfunc'),...
        fullfile(aap.directory_conventions.eeglabdir,'functions', 'resources'),...
        fullfile(aap.directory_conventions.eeglabdir,'functions', 'javachatfunc')...
        );
else
    % Check whether already in path, give warning if not
    if isempty(which('eeglab'))
       aas_log(aap,false,sprintf('EEG lab not found, if you need this you should add it to the matlab path manually, or set aap.directory_conventions.eeglabdir'));
    end;
end;

% Path to spm modifications to the top
addpath(fullfile(mfp,'extrafunctions','spm_mods'),'-begin');

%% Build required path list for cluster submission
% aa
reqpath=textscan(genpath(mfp),'%s','delimiter',':'); reqpath = reqpath{1};

p = textscan(path,'%s','delimiter',':'); p = p{1};

% spm
p_ind = cell_index(p,aap.directory_conventions.spmdir); % SPM-related dir
for ip = p_ind
    reqpath{end+1} = p{ip};
end
% spmtools
if isfield(aap.directory_conventions,'spmtoolsdir') && ~isempty(aap.directory_conventions.spmtoolsdir)
    SPMTools = textscan(aap.directory_conventions.spmtoolsdir,'%s','delimiter', ':');
    SPMTools = SPMTools{1};
    for pp = SPMTools'
        if exist(pp{1},'dir'), reqpath{end+1}=pp{1};end
    end
end

% MNE
if isfield(aap.directory_conventions,'mnedir') && ~isempty(aap.directory_conventions.mnedir)
    if exist(fullfile(aap.directory_conventions.mnedir,'matlab'),'dir')
        reqpath{end+1}=fullfile(aap.directory_conventions.mnedir,'matlab','toolbox');
        reqpath{end+1}=fullfile(aap.directory_conventions.mnedir,'matlab','examples');
    end
end

% EEGLAB
if ~isempty(aap.directory_conventions.eeglabdir)
    p_ind = cell_index(p,aap.directory_conventions.eeglabdir);
    for ip = p_ind
        reqpath{end+1} = p{ip};
    end
end;

% clean
reqpath=reqpath(strcmp('',reqpath)==0);
exc = cell_index(reqpath,'.git');
if exc, reqpath(exc) = []; end

aacache.reqpath = reqpath;
