% This script runs the Auditory fMRI example from the SPM manual
% (chapter 30, as of this writing) in aa. It uses BIDS version
% of the data, available for download at:
%
%  https://www.fil.ion.ucl.ac.uk/spm/download/data/MoAEpilot/MoAEpilot.bids.zip
%
% The dataset will be automatically downloaded if you have not done so yet.
function tutorial_2_SPM_CH30()

% variable names in ALLUPPERCASE are placeholders that
% you may want to edit before the script is run.
%
% make a copy of this file to edit and put it somewhere in your Matlab path

% -------------------------------------------------------------------------
% INITIALIZATION
% -------------------------------------------------------------------------
aa_ver5;

%% ------------------------------------------------------------------------
% Creating an aap structure
% -------------------------------------------------------------------------
% Call aarecipe to create the aap settings structure from a parameter file
% and a tasklist file.
% The default parameter file will be used (see tutorial_1_aa_setup), since
% no parameter filename is provided.
% The tasklist for tutorial_2_SPM_CH30 comes installed with aa.
%
% To use a non-default parameter file, set it as the first input to
% aarecipe, like
%PARAMETER_FNAME = '/path/to/parameter_xml_file';
%aap = aarecipe(PARAMETER_FNAME, tasklist_fname);

[aahome,~,~] = fileparts(which('aa_ver5'));
tasklist_fname = fullfile(aahome, 'examples', 'tutorials', 'tutorial_2_SPM_CH30.xml');
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
%
% When ROOT_PATH and .root setting are outcommented: use value from parameter file
%ROOT_PATH = '/path/to/dir/where/results_dir/will/be/created';
%aap.acq_details.root = ROOT_PATH;

% The analysisid is specific for a particular analysis, and typically set
% in the run script, i.e. here:
RESULTS_DIR = mfilename();
aap.directory_conventions.analysisid = RESULTS_DIR;

%% ------------------------------------------------------------------------
% Specifying the data directory
% -------------------------------------------------------------------------
% For BIDS data, point rawdatadir at the top level BIDS directory.
%
% In this tutorial, will use demo dataset MoAEpilot, that will be
% downloaded if necessary.
% Here it is assumed that aap.directory_conventions.rawdatadir in the
% parameter xml file is a single directory, not a list of directories
% separated by pathsep characters.
FULLDATAPATH = fullfile(aap.directory_conventions.rawdatadir, 'MoAEpilot');
aap.directory_conventions.rawdatadir = FULLDATAPATH;
aa_downloaddemo(aap, 'MoAEpilot');

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

end