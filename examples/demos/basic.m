% Automatic analysis (aa) - user master script
%
% This example runs one session of one subject through a standard SPM-based
% fMRI analysis. After running, look at the results using SPM in
% aa_demo_results/auditory/aamod_firstlevel_contrasts_00001/2014_03_29_9001/stats
%
% v3: Johan Carlin, MRC CBU, 2018-08-06
% v2: Tibor Auer MRC Cognition and Brain Sciences Unit, 2016-02-17
% v1: Rhodri Cusack Brain and Mind Institute, Western University, 2014-12-15

%% INIT
clear
aa_ver5;

%% LOAD TASKLIST
aap = aarecipe('basic.xml');

%% DEFINE STUDY SPECIFIC PARAMETERS
aap.options.wheretoprocess = 'localsingle'; % running locally

aap.tasksettings.aamod_norm_write.vox = [3 3 3];
aap.tasksettings.aamod_norm_write_meanepi.vox = [3 3 3];

%% DATA
% download the demo dataset (if necessary)
% Here it is assumed that aap.directory_conventions.rawdatadir in the
% parameter xml file is a single directory, not a list of directories
% separated by pathsep characters.
FULLDATAPATH = fullfile(aap.directory_conventions.rawdatadir, 'aa_demo');
aap.directory_conventions.rawdatadir = FULLDATAPATH;
aa_downloaddemo(aap, 'aa_demo');

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
aap.directory_conventions.analysisid = 'auditory';

% Add data
% Just one session
aap = aas_addsession(aap,'lullaby_task');
% Just one subject
aap = aas_addsubject(aap,'S1','2014_03_29_9001','functional',{6});

% Add model
% Just one regressor here ('Sound'): block onsets 0, 26, 52... 390 secs and duration 15 secs
aap = aas_addevent(aap,'aamod_firstlevel_model','*','*','Sound',0:26:390, 15);

% Specify contrast - just sound minus silence
aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts', '*', 'sameforallsessions', 1, 'sound-silence', 'T');

%% DO PROCESSING
aa_doprocessing(aap);
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
