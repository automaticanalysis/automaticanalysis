% Automatic analysis - initialise paths from recipe

function [aap]=aa_init(aap)

global aacache
aacache.bcp_path = path;

%% Set Paths
% Path for SPM
addpath(aap.directory_conventions.spmdir);
spm_jobman('initcfg');

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

% Path for EEGLAB
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

% Path to spm modifications to the top
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'extrafunctions','spm_mods'),'-begin');

%% Build required path list for cluster submission
% aa
mfp=textscan(mfilename('fullpath'),'%s','delimiter',filesep); mfp = mfp{1};
mfpi=find(strcmp('aa_engine',mfp));
reqpath=textscan(genpath([filesep fullfile(mfp{1:mfpi-1})]),'%s','delimiter',':'); reqpath = reqpath{1};

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
p_ind = cell_index(p,aap.directory_conventions.eeglabdir);
for ip = p_ind
    reqpath{end+1} = p{ip};
end

% clean
reqpath=reqpath(strcmp('',reqpath)==0);
reqpath(cell_index(reqpath,'.git')) = [];

aacache.reqpath = reqpath;
