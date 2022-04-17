# Example user scripts and task lists

The demos show more extended use cases. Additional requirements of the specific analyses in the demos and cbu directories are listed in the scripts.

#### basic.m & basic.xml
Runs one session of one subject through a standard SPM-based fMRI analysis.

#### branching1.m & branching1.xml
A variant of demo_basic demonstrating how branching can be used to explore how the order of slice time and motion correction affects results.

#### branching2.m & branching2.xml
A variant of demo_basic demonstrating nested multi-level branching.

#### ds000114.m & ds000114.xml
An example illustrating how to process a BIDS multimodal NIfTI dataset. Processing this dataset will benefit strongly from parallel processing with the 'batch' or 'parpool' queue processors.

#### ds000114_motor.m & ds000114_motor.xml
An simpler example analyzing only the motor task in ds000114

#### meeg.m & meeg.xml
Demonstrates a basic EEG pipeline on the LEMON dataset: http://fcon_1000.projects.nitrc.org/indi/retro/MPI_LEMON.html. Processing this dataset will benefit strongly from parallel processing with the 'batch' or 'parpool' queue processors.

#### denoising.m & denoising.xml
A demonstration of various motion correction strategies, including frame censoring, wavelet despiking, rWLS, and ICA
