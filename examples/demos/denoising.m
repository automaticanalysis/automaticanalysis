% a branched analysis to demonstrate the aa "denoising" modules, including:
%
%   1) FD- and DVARS-based scrubbing (frame censoring) w/ greyplot display
%   2) wavelet despiking, including epi rescaling
%   3) Jorn Diedrichsen's robust weighted least squares (rWLS)
%   4) FSL AROMA (i.e., untrained ICA denoising)
%
% data assumed: ds000114 motor task (10 subjects)
%
% notes
%
% 1) analysis is largely straightforward; some non-obvious stream renaming
% is required for AROMA (because the module re-runs some of the
% preprocessing using FSL equivalents) and rWLS. See stream
% customization code block below.
%
% 2) the rWLS and BrainWavelet toolboxes must be installed and toolbox
% entries included in the paramter file
%
% 2b)FSL and FSL/AROMA must be installed and the paths entered under the
% directory_conventions in the parameter file
%
%
%
% see comments in the tasklist (aaexample_denoising.xml) for
% (hopefully) helpful usage information


% ------------------------------------------------------------------------------------------------------------------------------
% INITIALIZATION
% ------------------------------------------------------------------------------------------------------------------------------

clear all;
aa_ver5;

aap = aarecipe('denoising.xml');      % this uses parameter file ~/.aa/aap_parameters_user.xml

% ------------------------------------------------------------------------------------------------------------------------------
% DIRECTORY AND DATA OPTIONS
% ------------------------------------------------------------------------------------------------------------------------------

ROOT_PATH = '/path/to/dir/where/results_dir/will/be/created';
RESULTS_DIR = 'name_of_results_directory';

aap.acq_details.root = ROOT_PATH;
aap.directory_conventions.analysisid = RESULTS_DIR;

% data specification
% for BIDS data, just point rawdatadir at the top level BIDS directory
% (i.e., wherever you downloaded ds000114)

FULLDATAPATH = '/full/path/to/toplevelBIDS';
aap.directory_conventions.rawdatadir = FULLDATAPATH;


aap.options.NIFTI4D = 1;
aap.options.autoidentifystructural_choosefirst = 1;
aap.options.autoidentifystructural_chooselast = 0;

% correct T1 effect using numdummies

aap.acq_details.numdummies = 4;
aap.acq_details.input.correctEVfordummies = 1;

% PCT customization - change this for your installation as appropriate
% (or just comment out to use default localsingle)

% aap.options.wheretoprocess='matlab_pct';
aap.options.wheretoprocess='parpool';   % better
aap.directory_conventions.poolprofile = 'local';
aap.options.aaparallel.numberofworkers = 15;

% ------------------------------------------------------------------------------------------------------------------------------
% STREAM CUSTOMIZATION
% ------------------------------------------------------------------------------------------------------------------------------

% rWLS docs recommend running on unsmoothed data, so redirect epi input
% (be sure to change 0003 to the correct ordinal number if you modify the tasklist...)

aap = aas_renamestream(aap,'aamod_firstlevel_model_00003','epi','aamod_norm_write_epi_00001.epi');

% explicitly feed coregged epi and structural to AROMA

aap = aas_renamestream(aap,'aamod_AROMA_denoise_00001','epi','aamod_coreg_extended_00001.epi');
aap = aas_renamestream(aap,'aamod_AROMA_denoise_00001','structural','aamod_coreg_extended_00001.structural');

% ------------------------------------------------------------------------------------------------------------------------------
% BIDS input and model
% ------------------------------------------------------------------------------------------------------------------------------

aap.acq_details.input.combinemultiple = true;

% you can use this to test a single subject (but secondlevel modeling will crash!)
% aap = aas_processBIDS(aap, [], {'finger_foot_lips'},{'sub-01'});

aap = aas_processBIDS(aap, [], {'finger_foot_lips'});

aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts_*', '*', 'sameforallsessions', [-0.5 -0.5 1], 'lips', 'T');

% ------------------------------------------------------------------------------------------------------------------------------
% RUN
% ------------------------------------------------------------------------------------------------------------------------------

aa_doprocessing(aap);

