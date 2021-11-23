
% This script runs the Auditory fMRI example from the SPM manual
% (chapter 30, as of this writing) in aa. It uses BIDS version
% of the data. available for download at:
%
%  https://www.fil.ion.ucl.ac.uk/spm/download/data/MoAEpilot/MoAEpilot.bids.zip
%
% (NB: this data is missing gzip'ed versions of the structral and epi which
% confuses aa. The easiest fix is just gzip /sub-01/anat/sub-01_T!w.nii
% and /sub-01/func/sub-01_task-auditory_bold.nii after downloading).

% make a copy of this file to edit and put it somewhere in your Matlab path

% -------------------------------------------------------------------------
% INITIALIZATION
% -------------------------------------------------------------------------

clear all;
aa_ver5;

% PARAMETER_FNAME is the aa parameter file you must customize
% before running aa -- see the aa documentation for help

PARAMETER_FNAME = '/path/to/parameter_xml_file';

% the tasklist comes installed with aa so you don't need to specify it

[aahome,~,~] = fileparts(which('aa_ver5'));
tasklist_fname = fullfile(aahome,'examples/tutorials/aa_example_SPM_CH30.xml');

aap = aarecipe(PARAMETER_FNAME,tasklist_fname);

% -------------------------------------------------------------------------
% results and data directory specification
% -------------------------------------------------------------------------

% results will be created in ROOT_PATH/RESULTS_DIR
% ROOT_PATH must exist on your machine
% aa will create RESULTS_DIR when it runs

ROOT_PATH = '/path/to/dir/where/results_dir/will/be/created';
RESULTS_DIR = 'name_of_results_directory';

aap.acq_details.root = ROOT_PATH;
aap.directory_conventions.analysisid = RESULTS_DIR;

% data specification
% for BIDS data, just point rawdatadir at the top level BIDS directory
% (i.e., wherever you downloaded MoAEpilot)

FULLDATAPATH = '/full/path/to/toplevelBIDS';

aap.directory_conventions.rawdatadir = FULLDATAPATH;

% -------------------------------------------------------------------------
% analysis options
% -------------------------------------------------------------------------

aap.options.autoidentifystructural_choosefirst = 1;
aap.options.autoidentifystructural_chooselast = 0;

aap.options.NIFTI4D = 1;
aap.acq_details.numdummies = 0;
aap.acq_details.input.correctEVfordummies = 0;

% make sure the report name is set
aap.directory_conventions.reportname='report.htm';


% -------------------------------------------------------------------------
% BIDS input
% -------------------------------------------------------------------------

aap = aas_processBIDS(aap);

% -------------------------------------------------------------------------
% modeling - contrast specification
% -------------------------------------------------------------------------

% UNITS can be 'secs' or 'scans' (the SPM auditory tutorial has it set
% for 'scans' in the manual but the BIDS tsv is in secs)

aap.tasksettings.aamod_firstlevel_model.xBF.UNITS = 'secs';

% processBIDS will create the events for the model, but you must define the contrasts

aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts', '*', 'sameforallsessions', 1, 'L_G_R','T');

% -------------------------------------------------------------------------
% run and report
% (you must have FSL installed to run reporting)
% -------------------------------------------------------------------------

aa_doprocessing(aap);
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
