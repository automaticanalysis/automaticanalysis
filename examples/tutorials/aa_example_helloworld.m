
% this is a simple aa script that runs a bare minimum tasklist using
% OpenNeuro dataset ds000114 (which you must download before running)
%
% This isn't a useful analysis; it's intended to check that your
% paths, parameter files, etc are set up correctly.
%
% See "aa_example_helloworld_model" to actually do something useful
% with the data.

% variable names in ALLUPPERCASE are placeholders that
% you must edit before the script can be run.

% make a copy of this file to edit and put it somewhere in your Matlab path

% -------------------------------------------------------------------------
% initialization
% -------------------------------------------------------------------------

clear all;
aa_ver5;

% PARAMETER_FNAME is the aa parameter file you must customize
% before running aa -- see the aa documentation for help

PARAMETER_FNAME = '/path/to/parameter_xml_file';

% the helloworld tasklist comes installed with aa:

[aahome,~,~] = fileparts(which('aa_ver5'));
tasklist_fname = fullfile(aahome,'examples/aa_example_helloworld.xml');

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
% BIDS input
% -------------------------------------------------------------------------

% we only run sub-01 for the finger_foot_lips task (ds000114
% contains four different tasks) since we're just testing whether
% your aa install is working...

aap = aas_processBIDS(aap, [], {'finger_foot_lips'}, {'sub-01'});

% -------------------------------------------------------------------------
% run
% -------------------------------------------------------------------------

aa_doprocessing(aap);
aa_close(aap);

% there really isn't any rea "results" to review other than
% checking that aa ran and created RESULTS_DIR and popluated 
% it with some files
