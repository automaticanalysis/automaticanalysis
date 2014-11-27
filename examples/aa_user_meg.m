clear;

%% Parameters

% Directories:
rawmri_dir    = '/megdata/camcan/camcan';
ana_dir    = '/imaging/ta02/aa/MEG';
ana_subdir = 'test';

%% AA PARAMETERS:

aa_ver4_nocloud

%% RECIPE:
aap = aarecipe('aap_parameters_defaults_CBSU.xml','aap_tasklist_meg.xml');

% SITE-SPECIFIC CONFIGURATION:
aap.options.wheretoprocess='qsub'; % queuing system			% typical value localsingle or qsub
aap.options.email='tibor.auer@mrc-cbu.cam.ac.uk';

%% DATA
% Directory & sub-directory for analysed data:
aap.acq_details.root = ana_dir;
aap.directory_conventions.analysisid = ana_subdir;

% Directory for raw data:
aap.directory_conventions.rawdatadir = rawmri_dir;

% Add subject (full):
aap = aas_add_meg_session(aap,'rest');
aap = aas_add_meg_session(aap,'task');
aap = aas_addsubject(aap,{[11 11] []},{'rest_raw.fif','task_raw.fif'});
aap = aas_addsubject(aap,{[11 78] []},{'rest_raw.fif','task_raw.fif'});

%% RUN
aa_doprocessing(aap);
aas_garbagecollection(aap,true);
