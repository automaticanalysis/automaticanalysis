# Example user scripts and task lists

The demos show more extended use cases. Additional requirements of the specific analyses in the demos and cbu directories are listed in the scripts.

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

#### aa_example_denoising
????