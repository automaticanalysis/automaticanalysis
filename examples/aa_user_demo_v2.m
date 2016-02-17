% Automatic analysis (aa) - user master script
% This example runs one session of one subject through a standard SPM-based
% fMRI analysis. After running, look at the results using SPM in
% ../analysis_v1/auditory/aamod_firstlevel_contrasts_00001/2014_03_29_9001/stats
% Prerequistes are aa and spm 8 or 12, I'm using Matlab 2014a.
% v2: Tibor Auer MRC Cognition and Brain Sciences Unit, 2016-02-17
% v1: Rhodri Cusack Brain and Mind Institute, Western University, 2014-12-15

%% INIT
clear

aa_ver4;

%% LOAD RECIPE AND TASKLIST
% You might want to set up a configuration for your local site settings in 
% aap_parameters_defaults_<your site>.xml and use it instead of aap_parameters_defaults.xml
% The settings for this script are manually configured in lines 20-42.
aap = aarecipe('aap_parameters_defaults.xml','aap_tasklist_demo.xml');
% SPM12 is used, set path according to your environment
aap.directory_conventions.spmdir = '/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r6685';
aap.directory_conventions.fsldir = '/imaging/local/software/fsl/v5.0.9/x86_64/fsl';
aap = aas_configforSPM12(aap);

%% DEFINE STUDY SPECIFIC PARAMETERS
aap.options.aa_minver = '4.3.0'; % designed for aa version 4.3.0 or above
aap.options.wheretoprocess = 'localsingle'; % running locally
aap.options.NIFTI4D = 1;

%% DATA
% Define following paths relative to this script. You might not always want
% to do this, but it is good for this example
localroot=fileparts(pwd);

% Location of raw DICOM data
aap.directory_conventions.rawdatadir = fullfile(localroot,'rawdata');
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
aap.acq_details.root = fullfile(localroot,'analysis_v2');
aap.directory_conventions.analysisid = 'auditory';     

% Add data
% Just one session
aap = aas_addsession(aap,'lullaby_task');
% Just one subject
aap = aas_addsubject(aap,'2014_03_29_9001',[6]);

% Add model
% Just one regressor here ('Sound'): block onsets 0, 26, 52... 390 secs and duration 15 secs 
aap = aas_addevent(aap,'aamod_firstlevel_model','*','*','Sound',0:26:390, 15);  

% Specify contrast - just sound minus silence
aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts', '*', 'sameforallsessions', 1, 'sound-silence', 'T');

%% DO PROCESSING
aap = aa_doprocessing(aap);
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));