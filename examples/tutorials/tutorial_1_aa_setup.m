% This is a script/function to help set up your user specific paramater xml
% file.
%
% It runs a bare minimum tasklist using OpenNeuro dataset ds000114 (that
% will be automatically downloaded if you have not done so yet).
% This isn't a useful analysis; it's intended to check that your
% paths, parameter files, etc are set up correctly.
%
% This file also shows the default function calls that will be used
% in an all aa run scripts / functions:
%   aa_ver5
%   aarecipe
%   aas_processBIDS % when using BIDS data
%   aa_doprocessing
%
% See "aa_example_helloworld_model" to actually do something useful
% with the data.
function tutorial_1_aa_setup()

% variable names in ALLUPPERCASE are placeholders that
% you may want to edit before the script is run.
%
% make a copy of this file to edit and put it somewhere in your Matlab path

%% ------------------------------------------------------------------------
% Initialization
% -------------------------------------------------------------------------
aa_ver5;

%% ------------------------------------------------------------------------
% Creating a parameter xml
% -------------------------------------------------------------------------
% Here we create a base parameter xml file, that will contain your default
% settings.
% It will be used each time aarecipe is called without explicitly passing a
% parameter filename.
% Settings can be overriden later where needed, in specific parameter files
% and in run scripts, but that is not a subject in this tutorial file.
if ~aa_has_user_parameter_file()
    aas_create_parameter_xml('', true, 'use_default_location', true, 'use_default_filename', true);
end

%% ------------------------------------------------------------------------
% Creating an aap structure
% -------------------------------------------------------------------------
% Call aarecipe to createa the aap settings structure from a parameter file
% and a tasklist file.
% The default parameter file created in the previous step will be used,
% since no parameter filename is provided.
% The helloworld tasklist comes installed with aa:
[aahome,~,~] = fileparts(which('aa_ver5'));
tasklist_fname = fullfile(aahome,'examples' ,'tutorials', 'tutorial_1_aa_setup.xml');
aap = aarecipe(tasklist_fname);

%% ------------------------------------------------------------------------
% Specifying the results directories
% -------------------------------------------------------------------------
% Results will be created in the directory specified in
% aap.acq_details.root, in a specific subdirectory as specified in
% aap.directory_conventions.analysisid
% It is possible within a run script to override the default values taken
% from the parameter file.
%
% The directory in aap.acq_details.root must exist on your machine
% When it runs, aa will create a subdirectory RESULTS_DIR.
%aap.acq_details.root = ROOT_PATH; % Outcommented: use value from parameter file
RESULTS_DIR = mfilename(); % Typically specific for a particular analysis.
aap.directory_conventions.analysisid = RESULTS_DIR;

%% ------------------------------------------------------------------------
% Specifying the data directory
% -------------------------------------------------------------------------
% For BIDS data, point rawdatadir at the top level BIDS directory.
%
% In this tutorial, will use a demo dataset ds000114, that will be
% downloaded if necessary.
% Here it is assumed that aap.directory_conventions.rawdatadir in the
% parameter xml file is a single directory, not a list of directories
% separated by pathsep characters.
FULLDATAPATH = fullfile(aap.directory_conventions.rawdatadir, 'ds000114');
aap.directory_conventions.rawdatadir = FULLDATAPATH;
aa_downloaddemo(aap, 'ds000114');

%% ------------------------------------------------------------------------
% Using BIDS input
% -------------------------------------------------------------------------
% When using BIDS input, run aas_processBIDS before aa_doprocessing
%
% Here we only run sub-01 for the finger_foot_lips task (ds000114
% contains four different tasks) since we're just testing whether
% your aa install is working.
aap = aas_processBIDS(aap, [], {'finger_foot_lips'}, {'sub-01'});

%% ------------------------------------------------------------------------
% Run
% -------------------------------------------------------------------------
aa_doprocessing(aap);

% There really isn't any real "results" to review other than
% checking that aa ran and created RESULTS_DIR and populated
% it with some files.

end