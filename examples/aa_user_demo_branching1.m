% Automatic analysis (aa) - user master script This example runs one session of one
% subject through a standard SPM-based fMRI analysis. See also aa_user_demo.
%
% This script demonstrates how branching can be used to explore how the order of slice
% time and motion correction affects results (see aap_tasklist_demo_branching1.xml).
%
% v3: Johan Carlin, MRC CBU, 2018-08-06
% v2: Tibor Auer MRC Cognition and Brain Sciences Unit, 2016-02-17
% v1: Rhodri Cusack Brain and Mind Institute, Western University, 2014-12-15

%% INIT
clear
aa_ver5;

%% LOAD TASKLIST
aap = aarecipe('aap_tasklist_demo_branching1.xml');

%% DEFINE STUDY SPECIFIC PARAMETERS
aap.options.wheretoprocess = 'localsingle';

aap.tasksettings.aamod_slicetiming(1).sliceorder = [1:2:36 2:2:36];
aap.tasksettings.aamod_slicetiming(1).refslice = 16;
aap.tasksettings.aamod_norm_write(1).vox = [3 3 3];
aap.tasksettings.aamod_norm_write_meanepi(1).vox = [3 3 3];
aap.tasksettings.aamod_slicetiming(2).sliceorder = [1:2:36 2:2:36];
aap.tasksettings.aamod_slicetiming(2).refslice = 16;
aap.tasksettings.aamod_norm_write(2).vox = [3 3 3];
aap.tasksettings.aamod_norm_write_meanepi(2).vox = [3 3 3];

% download the demo dataset (if necessary)
aap = aa_downloaddemo(aap);

% Define how subject identifier (e.g. 2014_03_29_9001) is turned into
% subject foldername in rawdatadir
aap.directory_conventions.subjectoutputformat = '%s';
aap.directory_conventions.dicomfilter='*.IMA';
% This is name of the structural protocol, as typed into the scanner. If 
% you're consistent with this, it can be found automatically in each new subject
aap.directory_conventions.protocol_structural = 'MPRAGE  iPAT2_sag';

% Number of dummy scans at start of EPI runs
aap.acq_details.numdummies = 10;

%% STUDY
% Where to put the analyzed data
aap.acq_details.root = fullfile(aap.acq_details.root,'aa_demo');
aap.directory_conventions.analysisid = 'auditory_branching1';     

% Add data
% Just one session
aap = aas_addsession(aap,'lullaby_task');
% Just one subject
aap = aas_addsubject(aap,'S1','2014_03_29_9001','functional',{6});

% Add model
% Just one regressor here ('Sound'): block onsets 0, 26, 52... 390 secs and duration 15 secs 
aap = aas_addevent(aap,'aamod_firstlevel_model_*','*','*','Sound',0:26:390, 15);  

% Specify contrast - just sound minus silence
aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts_*', '*', 'sameforallsessions', 1, 'sound-silence', 'T');

%% DO PROCESSING
aa_doprocessing(aap);
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
