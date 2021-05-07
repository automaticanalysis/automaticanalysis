# Example user scripts and task lists
This directory contains example scripts that use publicly available data (either BIDS or Cusack lab). Requirements of the specific analyses beyond those listed on the main project page are listed in the scripts. For an automated test, see aa_testcode/aatest. 

#### demo_basic_aa.m & demo_basic_tasklist.xml
Runs one session of one subject through a standard SPM-based fMRI analysis.

#### demo_branching1_aa.m & demo_branching1_tasklist.xml
A variant of demo_basic demonstrating how branching can be used to explore how the order of slice time and motion correction affects results.

#### demo_branching2_aa.m & demo_branching2_tasklist.xml
A variant of demo_basic demonstrating nested multi-level branching.

#### bids_ds000114_aa.m & bids_ds000114_tasklist.xml
An example illustrating how to process a BIDS multimodal NIfTI dataset. Processing this dataset will benefit strongly from parallel processing with the 'batch' or 'parpool' queue processors.

#### meeg_aa.m & meeg_tasklist.xml
Demonstrates a basic EEG pipeline on the LEMON dataset: http://fcon_1000.projects.nitrc.org/indi/retro/MPI_LEMON.html. Processing this dataset will benefit strongly from parallel processing with the 'batch' or 'parpool' queue processors.

#### aa_downloaddemo.m
Downloads dataset from user-specified URL into aa_demo folder as specified in  aap.directory_conventions.rawdatadir. If no URL is provided then the small aa_demo from Rhodri Cusack will be used.

#### aamod_template_session (.m & .xml)
 Templates for code and header files for new MRI analysis modules.
