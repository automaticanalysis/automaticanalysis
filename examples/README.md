# Example user scripts and task lists
This directory contains example scripts that use publicly available data (either BIDS or Cusack lab). For an automated test, see aa_testcode/aatest.

Besides the requirements listed on the main project page, the following extensions/materials should be present/installed on your system:

EEGLab plugins:
- "AMICA"
- "Fieldtrip-lite"
- "Fileio"
- "ICLabel"
- "Viewprops"
- "clean_rawdata"
- "dipfit"
- "firfilt"
- "fitTwoDipoles"
- "neuroscanio"

Specifically for the meeg example:
Please read the instructions on the template directory you as the user have to set up in  https://github.com/fieldtrip/fieldtrip/blob/master/bin/ft_postfreesurferscript.sh and assemble all files accordingly. The path to the resulting template directory must be known to aa (e.g. listed in aa_parameters_user).

#### aa_user_demo (.m & .xml)
Runs one session of one subject through a standard SPM-based fMRI analysis.

#### aa_user_demo_branching1 (.m & .xml)
A variant of aa_user_demo demonstrating how branching can be used to explore how the order of slice time and motion correction affects results.

#### aa_user_demo_branching2 (.m & .xml)
A variant of aa_user_demo demonstrating  nested multi-level branching.

#### aa_user_bids_ds000114 (.m & .xml)
An example illustrating how to process a BIDS multimodal NIfTI dataset. Processing this dataset will benefit strongly from parallel processing with the 'batch' or 'parpool' queue processors.

#### aa_user_meeg (.m & .xml)
Demonstrates a basic EEG pipeline on the LEMON dataset: http://fcon_1000.projects.nitrc.org/indi/retro/MPI_LEMON.html. Processing this dataset will benefit strongly from parallel processing with the 'batch' or 'parpool' queue processors.

#### aa_downloaddemo.m
 Downloads dataset from user-specified URL into aa_demo folder as specified in  aap.directory_conventions.rawdatadir. If no URL is provided then the small aa_demo from Rhodri Cusack will be used.

#### aamod_template_session (.m & .xml)
 Templates for code and header files for new MRI analysis modules.