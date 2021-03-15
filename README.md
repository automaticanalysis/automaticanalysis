Automatic Analysis software was originally developed by Dr Rhodri Cusack
for supporting research at the MRC Cognition and Brain Science Unit. It
is made available to the academic community in the hope that it may
prove useful.

[![Docker image available!](https://img.shields.io/badge/docker-bids%2Faa-brightgreen.svg?logo=docker&style=flat)](https://hub.docker.com/r/bids/aa/)
[![manuscript](https://img.shields.io/badge/doi-10.3389%2Ffninf.2014.00090-blue.svg?style=flat)](http://doi.org/10.3389/fninf.2014.00090)
[![resource](https://img.shields.io/badge/RRID-SCR__003560-blue?style=flat)](https://scicrunch.org/scicrunch/Resources/record/nlx_144509-1/SCR_003560/resolver)

Definitions: aa means the Automatic Analysis software package and any
associated documentation whether electronic or printed.

# License

Use of this software is subject to the terms of the license, found in
the license.txt file distributed with this software.

# Update #

Please, be aware that this README is not updated anymore; therefore, some information might be outdated and the links might be broken.

The new website at http://automaticanalysis.github.io is now live and should be the primary point of reference.

Please, consult with this README only when you do not find the desired information at http://automaticanalysis.github.io.

# About

aa is a pipeline system for neuroimaging written primarily in Matlab. It
supports SPM 12, as well as selected functions from other software packages 
(e.g., FSL, Freesurfer). The goal is to facilitate automatic, flexible, and 
replicable neuroimaging analyses through a comprehensive pipeline system.

More information can be found on the [wiki](https://github.com/automaticanalysis/automaticanalysis/wiki).

# Requirements

 - System: Linux or MacOX (Windows is not supported. Sorry!)
 - Softwares: 
   - **MATLAB** - It has been tested with version r2013a and later
   - **SPM** - It has been tested with versions SPM12 and later (SPM12 r7487 or later **required**!)
 - Recommended software for features:
   - For MEG (maxfilter): Neuromag (Elekta Instrument AB Stockholm, Stockholm, Sweden)
   - For visual representation of the pipeline: [GraphViz](http://www.graphviz.org)
 - Recommended/Partially supported software for analysis:
   - [Advanced Normalization Tools](http://stnava.github.io/ANTs)
   - [VBM8](http://www.neuro.uni-jena.de/vbm)
   - [FMRIB Software Library](http://fsl.fmrib.ox.ac.uk/fsl/fslwiki)
   - [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/fswiki)
   - [Automatic Defacing Tools templates (*.gca)](https://surfer.nmr.mgh.harvard.edu/fswiki/mri_deface)
   - [Group ICA Of fMRI Toolbox](http://mialab.mrn.org/software/gift/index.html)
   - [EEGLab](http://sccn.ucsd.edu/eeglab)
   - [BrainWavelet Toolbox](http://www.brainwavelet.org)
   - [Motion FingerPrint](https://www.medizin.uni-tuebingen.de/kinder/en/research/neuroimaging/software)
 
# Download

aa is maintained on [github](https://github.com/automaticanalysis/automaticanalysis/).

# Install

Unzip the downloaded package, and add the main folder to the MATLAB path.

# Configure

Some parameters (path settings, format settings) stored in **_<aarootdir>/aa_parametersets/aap_parameters_defaults.xml_** have to be configured to reflect your local environment. We recommend naming your local configuration file aap_parameters_user.xml if you want to run the examples easily (see below).
 
 Instead of copying the main **_aap_parameters_defaults.xml_**, we recommend to create a local version by *inheriting* most of the settings from **_aap_parameters_defaults.xml_**. See **_aap_parameters_defaults_CBSU.xml_** as an example.
 - *directory_conventions/rawdatadir*: Directories to find raw MRI data
 - *directory_conventions/subjectoutputformat*: `sprintf` formatting string to get subject directory as stored in *directory_conventions/rawdatadir*
 - *directory_conventions/seriesoutputformat*: `sprintf` formatting string to get series directory as stored in subject directory
 - *directory_conventions/protocol_structural*: Name of the structural/anatomical protocol as stored in the DICOM header
 - *directory_conventions/dicomfilter*: Directory listing filter to find DICOM data
 - *directory_conventions/toolboxes/spm*: Settings for SPM. **N.B.: You should not modify SPM version in your user script but rather in your local parameter file.**
Optional:
 - For distortion correction using fieldmap:
   - *directory_conventions/protocol_fieldmap*: Name of the fieldmap protocol as stored in the DICOM header

 - For multichannel segmentation using T2-weighted image:
   - *directory_conventions/protocol_t2*: Name of the T2-weighted protocol as stored in the DICOM header

 - For MEEG: 
   - *directory_conventions/rawmeegdatadir*: Directories to find raw MEEG data
   - *directory_conventions/meegsubjectoutputformat*: `sprintf` formatting string to get subject directory as stored in *directory_conventions/rawmeegdatadir*		
   - (optional) *directory_conventions/neuromagdir*: Path to Neuromag (maxfilter)

 - Optional analysis software:
   - FMRIB Software Library:
     - *directory_conventions/fsldir*: Path to FSL
     - *directory_conventions/fslsetup*: Path to FSL setup/configuration script to be executed before any FSL command
     - *directory_conventions/fslshell*: Shell used to run FSL. It should be in accordance with *directory_conventions/fslsetup*
     - *directory_conventions/fsloutputtype*: Output type used by FSL ("NIFTI" recommended)

   - FreeSurfer:
     - *directory_conventions/freesurferdir*: Path to FreeSurfer
     - *directory_conventions/freesurfersetup*: Path to FreeSurfer setup/configuration script, executing before any FreeSurfer command
     - *directory_conventions/freesurfershell*: Shell used to run FreeSurfer. It should be in accordance with *directory_conventions/freesurfersetup* and *directory_conventions/freesurferenvironment*
     - *directory_conventions/freesurferenvironment*: Path to FreeSurfer environmental setup script, executing before any FreeSurfer command
   
   - *directory_conventions/toolboxes/eeglab*: Settings for EEGLab
   - *directory_conventions/GIFTdir*: Path to GIFT
   - *directory_conventions/ANTSdir*: Path to ANTS
   - *directory_conventions/toolboxes/fieldtrip*: Settings for FieldTrip
   - *directory_conventions/BrainWaveletdir*: Path to BrainWavelet
   - *directory_conventions/spmtoolsdir*: Path to external/custom SPM toolboxes

 - Parallel processing: 
   - *directory_conventions/poolprofile*: MATLAB function to create cluster profile for parallel processing in a local job (should be in your MATLAB path)

# Get started
The user scripts and task lists in the examples folder demonstrate a wide range of analyses. It is easy to get started with AA by adapting the example that most resembles your planned analysis. Before you do that, we recommend that you try to run a few these scripts to confirm that your installation is running (e.g., aa_user_demo.m).

# Testing
We provide rudimentary automated tests under the aa_test folder. If you are having trouble with AA it's a good idea to start by checking that you can get the automated test to run.



# Help and support

Please feel free to open an issue here on Github if you are having trouble getting AA to work.

There is also a discussion list:

https://groups.google.com/d/forum/automaticanalysis

# Software updates

The master branch on github is only updated after testing, so this version is the one to use. Developer versions are maintained in each developer's AA fork on Github.

# References and citation

For any papers that report data analyzed with aa, please include the
github repo URL and cite the aa paper:

Cusack R, Vicente-Grabovetsky A, Mitchell DJ, Wild CJ, Auer T, Linke AC,
Peelle JE (2015) Automatic analysis (aa): Efficient neuroimaging
workflows and parallel processing using Matlab and XML. Frontiers in
Neuroinformatics 8:90.
http://dx.doi.org/10.3389/fninf.2014.00090
