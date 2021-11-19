
% this script runs basic preprocessing, and first- and second-level model
% on the ds000114 dataset used in the "helloworld" example script. It also
% demonstrates how to do processing using the Matlab Parallel Computing
% Toolbox if it is available.

% variable names in ALLUPPERCASE are placeholders that
% you must edit before the script can be run

% make a copy of this file to edit and put it somewhere in your Matlab path

% -------------------------------------------------------------------------
% INITIALIZATION
% -------------------------------------------------------------------------

more off;
clear all;
aa_ver5;

% PARAMETER_FNAME is the aa parameter file you must customize
% before running aa -- see the aa documentation for help

PARAMETER_FNAME = '/path/to/parameter_xml_file';

% the helloworld tasklist comes installed with aa:

[aahome,~,~] = fileparts(which('aa_ver5'));
tasklist_fname = fullfile(aahome,'examples/tutorials/aa_example_helloworld_model.xml');

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
% (i.e., wherever you downloaded ds000114)

FULLDATAPATH = '/full/path/to/toplevelBIDS';

aap.directory_conventions.rawdatadir = FULLDATAPATH;

% -------------------------------------------------------------------------
% analysis options
% -------------------------------------------------------------------------

% there are many (many) aa analysis options
% here are just a few examples

% there are multiple structuals in the ds000114 dataset
% we need to tell aa to use the first one

aap.options.NIFTI4D = 1;
aap.options.autoidentifystructural_choosefirst = 1;
aap.options.autoidentifystructural_chooselast = 0;

% correct for T1 effect using numdummies

aap.acq_details.numdummies = 4;
aap.acq_details.input.correctEVfordummies = 1;

% make sure the report name is set
aap.directory_conventions.reportname='report.htm';


% >>>>>>>> Parallel Computing Toolbox <<<<<<<<<
% you will need to optimize these settings for your processor

% to test whether the PCT is available on your machine, check if you have
% the Parallel Computing Toolbox installed using "ver" If so, check
% feature('numcores') -- if you don't have more than, say, 4 cores,
% then using the PCT isn't going to help you.

% aap.options.wheretoprocess='matlab_pct';
% aap.directory_conventions.poolprofile = 'local';
% aap.options.aaparallel.numberofworkers = 15;

% -------------------------------------------------------------------------
% BIDS input
% -------------------------------------------------------------------------

% ds000114 has the data structured as two separate sessions which
% we want to combine into a single analysis. (This is somewhat atypical
% in BIDS data.) We do so using the combinemultiple option. Note aa may
% print a warning about this -- you can ignore it.

aap.acq_details.input.combinemultiple = true;

% Here are a three options for processing the data:

% 1) this will run all subjects and all tasks. Since the contrast
% specified in the next section only works for the ds000114 motor
% task (finger_foot_lips), we don't do that:

% aap = aas_processBIDS(aap);

% 2) this will run all subjects for the specified task:

aap = aas_processBIDS(aap, [], {'finger_foot_lips'});

% 3) you can run only one subject to test the firstlevel modeling,
% but aamod_secondlevel_model will likely crash without at
% least four subjects:

%  aap = aas_processBIDS(aap, [], {'finger_foot_lips'}, {'sub-01'});

% -------------------------------------------------------------------------
% modeling
% -------------------------------------------------------------------------

% aa will read the event files in ds000114 and create a model, but BIDS
% currently does not support contrast specification. As such model
% contrasts must be added here.

% note aas_addcontrast must come AFTER aas_processBIDS in the mfile

aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts_*', '*', 'sameforallsessions', [-0.5 -0.5 1], 'lips', 'T');
aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts_*', '*', 'singlesession:finger_foot_lips_test', [-0.5 -0.5 1], 'lips-test', 'T');
aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts_*', '*', 'singlesession:finger_foot_lips_retest', [-0.5 -0.5 1], 'lips-retest', 'T');

% -------------------------------------------------------------------------
% run
% -------------------------------------------------------------------------

aa_doprocessing(aap);

% aa_report will crawl the results and generate an HTML summary
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));

aa_close(aap);

