% FILE: SPM_CH30.m
%
% This script runs the Auditory fMRI example from the SPM manual
% (chapter 30, as of this writing) in aa. 
%
% Make a copy of this file to edit and put it somewhere in your Matlab path
% along with a copy of the tasklist SPM_CH30.xml
%
% Variable names in ALLUPPERCASE are placeholders that you will need to
% customize before the script can be run.

% Data
%
% This script uses the BIDS version of the ch 30 data, available at:
%
%  https://www.fil.ion.ucl.ac.uk/spm/download/data/MoAEpilot/MoAEpilot.bids.zip
%
% Download & unzip this file and place the folder somewhere convenient
% on your machine. Alternatively, aa  will attempt to automatically
% download the data for you when the script runs (see comment at
% aa_downloaddemo below).

% -------------------------------------------------------------------------
% 0) aa initialization
% -------------------------------------------------------------------------

% a call to aa_ver5 is required as the first line of an aa script

aa_ver5;

% -------------------------------------------------------------------------
% 1) initializing the aap structure
% -------------------------------------------------------------------------
%
% An aa script begins with a call to aarecipe to create an aap "settings" 
% structure. This function takes a tasklist file and (optionally) a
% parameter file. If no parameter file is passed, the default parameter 
% file will be used.
%
% The default parameter file is assumed to be:
%
%       $HOME/.aa/aap_parameters_user.xml
%
% This file must exist if a parameter file is not passed to aarecipe.
%
% To use a named parameter file (customized, say, for a specific analysis)
% pass a path to the file as the first parameter in aarecipe:
%
%   PARAMETER_FNAME = '/path/to/parameter_file';
%   aap = aarecipe(PARAMETER_FNAME, tasklist_fname);
%

% Note we assued the tasklist is named "SPM_CH30.xml". It is common
% practice to keep this file in the same directory as the "userscript"
% (the file you are currently viewing) and having the the same name 
% (albeit with an xml extension).

aap = aarecipe('SPM_CH30.xml');

% ------------------------------------------------------------------------
% 2) specify the data directory
% -------------------------------------------------------------------------
%
% aa will look for data to be used in a given analysis in 
%
%   aap.directory_conventions.rawdatadir
%
% when using BIDS data, this is simply the top level BIDS directory.
%
% If you downloaded the MoAEpilot data, set DATAPATH to the fullpath
% name of where the data is located:

% for example: DATA_PATH = '/volumes/bigdisk/imaging/MoAEpilot'

DATA_PATH = '/fullpath/to/MoAEpilot';
aap.directory_conventions.rawdatadir = DATA_PATH;

% If you would like to have aa attempt to automatically download the 
% data for you, set autodownloadflag to true

autodownloadflag = false;

if (autodownloadflag == true)
    aa_downloaddemo(aap, 'MoAEpilot');
end

% note that automatic data download occasionally fails due to
% network issues, server availablilty, etc)

% ------------------------------------------------------------------------
% 3) specify the results directory
% -------------------------------------------------------------------------
%
% aa will save analysis results to the directory:
%
%   aap.acq_details.root/aap.directory_conventions.analysisid
%
% aap.acq_details.root must exist (you must create it) but aa
% will create aap.directory_conventions.analysisid if need be
%
% Some aa users prefer to always use the default settings of these 
% parameters (set in the parameter file). If you used the aa parameter file
% utility, you specified these directories during setup). However, here we
% will explicity set these parameters for illustrative purposes (parameters
% set in the userscript override the value set in the parameter file).

% for example: ROOT_PATH = '/volumes/bigdisk/imaging'
%              RESULTS_dir = 'MoAEpilot_RESULTS'

% (it is generally good practice to keep results separate from the data)

ROOT_PATH = '/path/to/dir/where/results_directory/will/be/created';
aap.acq_details.root = ROOT_PATH;

RESULTS_DIR = 'name_of_results_directory';
aap.directory_conventions.analysisid = RESULTS_DIR;

% -------------------------------------------------------------------------
% 4) specify analysis options
% -------------------------------------------------------------------------

% here we demonstrate the bare minimum of customizing analysis settings

% note variables not set here will use default values specified either 
% in your parameter file or taken from the module header (see any file
% aamod_*.xml in the "aa_modules" directory aa distribution for examples)

% the BIDS version of the auditory tutorial data is a 4D nifti file:

aap.options.NIFTI4D = 1;

% the SPM manual describes removing inital volumes because of T1
% effects. Howver, these volumes have been omitted when the BIDS
% data was generated so adjustment is unnecessary here. Otherwise we would 
% set numdummies > 0 and also set correctEVfordummies = 1 (true) so that
% event timing (read from the BIDS tsv files) would be corrected.

aap.acq_details.numdummies = 0;
aap.acq_details.input.correctEVfordummies = 0;

% the SPM manual specifies a smoothing kernal of 6 mm. We can set this
% here in the aap struct. However, for illustrative purposes we set
% this parameter in the tasklist (see aamod_smooth in SPM_CH30.xml).
% Either approach is valid.

% aap.tasksettings.aamod_smooth.FWHM = 6;

% UNITS can be 'secs' or 'scans' (the SPM auditory tutorial has it set
% for 'scans' in the manual but a BIDS tsv is always specified in secs)

aap.tasksettings.aamod_firstlevel_model.xBF.UNITS = 'secs';

% pick a name for the analysis report (this file is written in the
% results directory when reporting is complete)

aap.directory_conventions.reportname='report.htm';

% -------------------------------------------------------------------------
% 5) process BIDS input
% -------------------------------------------------------------------------

% aas_processBIDS inputs and processes the BIDS data

aap = aas_processBIDS(aap);

% -------------------------------------------------------------------------
% 6) modeling - contrast specification
% -------------------------------------------------------------------------

% aas_processBIDS will define the events for the model (these are read from
% the tsv files in the BIDS directory), but you must define the contrasts 
% that appear in your model using aas_addcontrast

% note any calls to aas_addcontrast must appear *after* aas_processBIDS

aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts', '*', 'sameforallsessions', 1, 'L_G_R','T');

% -------------------------------------------------------------------------
% 7) run and report
% -------------------------------------------------------------------------

aa_doprocessing(aap);
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));

% done!
