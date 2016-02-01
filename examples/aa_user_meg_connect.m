% This is an example how to connect MEG
% It connects to a fully preprocessed dataset in /imaging/xy01/aa/BP
% It uses a tasklist aap_tasklist_meg_epochs.xml (not included) with a single main module aamod_meg_epochs
%
% Tibor Auer, MRC-CBSU
% 01-02-2016

clear
aa_ver4

%% RECIPE:
aap = aarecipe('aap_parameters_defaults_CBSU.xml','aap_tasklist_meg_epochs.xml');

% SITE-SPECIFIC CONFIGURATION:
aap.options.wheretoprocess='localsingle'; % queuing system			% typical value localsingle or qsub
aap.options.email='xy01@mrc-cbu.cam.ac.uk';
aap.tasksettings.aamod_meg_epochs.timewindow = [-2000 500];

%% DATA
% Directory & sub-directory for analysed data:
aap.acq_details.root = '/imaging/xy01/aa';
aap.directory_conventions.analysisid = 'MEGconnect';

% Directory for raw data:
aap.directory_conventions.rawmegdatadir = '/megdata/cbu/ftd/';

% Add subject (full):
aap = aas_add_meg_session(aap,'psp_bp');
aap = aas_addsubject(aap,{[12 442] []},{'psp_button_press_self_raw.fif'});

% Add conditions
aap = aas_add_meg_event(aap,'aamod_meg_epochs',[12 442],'psp_bp','BP',8192,34);

% Connect
aap=aas_doprocessing_initialisationmodules(aap);
aap.directory_conventions.allowremotecache = 0;
remotePipes = struct('host',           '', ...
    'directory',      '/imaging/xy01/aa/BP', ...
    'allowcache',     0, ...
    'maxstagetag',   'aamod_meg_denoise_ICA_2_applytrajectory_00001', ...
    'checkMD5',       1);
aap = aas_connectAApipelines(aap, remotePipes);

%% RUN
aa_doprocessing(aap);