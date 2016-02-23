Automated Analysis software was originally developed by Dr Rhodri Cusack
for supporting research at the MRC Cognition and Brain Science Unit. It
is made available to the academic community in the hope that it may
prove useful.

Definitions: aa means the Automated Analysis software package and any
associated documentation whether electronic or printed.

# License

Use of this software is subject to the terms of the license, found in
the license.txt file distributed with this software.

# About

aa is a pipeline system for neuroimaging written primarily in Matlab. It
robustly supports recent versions of SPM, as well as selected functions
from other software packages. The goal is to facilitate automatic,
flexible, and replicable neuroimaging analyses through a comprehensive
pipeline system.

More information can be found on the aa webpage:

http://www.automaticanalysis.org

# Requirements

 - System: Linux or MacOX (Windows is not supported. Sorry!)
 - Softwares: 
   - **MATLAB** - It has been tested with version r2009a and later
   - **SPM** - It has been tested with versions SPM5 and later (SPM12 is recommended!)
 - Recommended software for features:
   - For MEG (maxfilter): Neuromag (Elekta Instrument AB Stockholm, Stockholm, Sweden)
   - For visual representation of the pipeline: [GraphViz](http://www.graphviz.org)
   - For generating and serialising provenance: [Raptor RDF Syntax Library](http://librdf.org/raptor)
 - Recommended/Partially supported software for analysis:
   - [Advanced Normalization Tools](http://stnava.github.io/ANTs)
   - [VBM8](http://www.neuro.uni-jena.de/vbm)
   - [FMRIB Software Library](http://fsl.fmrib.ox.ac.uk/fsl/fslwiki)
   - [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/fswiki)
   - [Automated Defacing Tools templates (*.gca)](https://surfer.nmr.mgh.harvard.edu/fswiki/mri_deface)
   - [Group ICA Of fMRI Toolbox](http://mialab.mrn.org/software/gift/index.html)
   - [EEGLab](http://sccn.ucsd.edu/eeglab)
   - [Advanced Normalization Tools](http://picsl.upenn.edu/software/ants)
   - [BrainWavelet Toolbox](http://www.brainwavelet.org)
   - [Motion FingerPrint](https://www.medizin.uni-tuebingen.de/kinder/en/research/neuroimaging/software)
 
# Download

aa is maintained on github. The most recent version can be found at:

https://github.com/rhodricusack/automaticanalysis/

# Install

Unzip the downloaded package, and add the main folder to the MATLAB path.

# Configure

Some parameters (path settings, format settings) stored in **_<aarootdir>/aa_recipes_and_parametersets/aap_parameters_defaults.xml_** have to be configured. Instead of editing the main **_aap_parameters_defaults.xml_**, we recommend to create a local version by *inhereting* most of the settings from **_aap_parameters_defaults.xml_**. See **_aap_parameters_defaults_CBSU.xml_** as an example.
 - *directory_conventions/rawdatadir*: Directories to find raw MRI data
 - *directory_conventions/subjectoutputformat*: `sprintf` formatting string to get subject directory as stored in *directory_conventions/rawdatadir*
 - *directory_conventions/seriesoutputformat*: `sprintf` formatting string to get series directory as stored in subject directory
 - *directory_conventions/protocol_structural*: Name of the structural/anatomical protocol as stored in the DICOM header
 - *directory_conventions/dicomfilter*: Directory listing filter to find DICOM data
 - *directory_conventions/spmdir*: Path to SPM
Optional:
 - For distortion correction using fieldmap:
   - *directory_conventions/protocol_fieldmap*: Name of the fieldmap protocol as stored in the DICOM header

 - For multichannel segmentation using T2-weighted image:
   - *directory_conventions/protocol_t2*: Name of the T2-weighted protocol as stored in the DICOM header

 - For MEG: 
   - *directory_conventions/rawmegdatadir*: Directories to find raw MEG data
   - *directory_conventions/megsubjectoutputformat*: `sprintf` formatting string to get subject directory as stored in *directory_conventions/rawmegdatadir*		
   - *directory_conventions/neuromagdir*: Path to Neuromag (maxfilter)

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
   
   - *directory_conventions/eeglabdir*: Path to EEGLab
   - *directory_conventions/GIFTdir*: Path to GIFT
   - *directory_conventions/ANTSdir*: Path to ANTS
   - *directory_conventions/fieldtripdir*: Path to FieldTrip
   - *directory_conventions/BrainWaveletdir*: Path to BrainWavelet
   - *directory_conventions/spmtoolsdir*: Path to external/custom SPM toolboxes

 - Parallel processing: 
   - *directory_conventions/qsubscheduler*: MATLAB function to set up qsub scheduler to execute jobs in the cluster (should be in your MATLAB path)
   - *directory_conventions/poolprofile*: MATLAB function to create cluster profile for parallel processing in a local job (should be in your MATLAB path)

# Help and support

In addition to the wiki, there is also a discussion list:

https://groups.google.com/d/forum/automaticanalysis

# Software updates

The master branch on github will always have the latest stable version.
To be notified when a new version is available, please sign up for the
automaticanalysis-announce email list:

https://groups.google.com/d/forum/automaticanalysis-announce

# References and citation

For any papers that report data analyzed with aa, please include the
website (http://www.automaticanalysis.org) and cite the aa paper:

Cusack R, Vicente-Grabovetsky A, Mitchell DJ, Wild CJ, Auer T, Linke AC,
Peelle JE (2015) Automatic analysis (aa): Efficient neuroimaging
workflows and parallel processing using Matlab and XML. Frontiers in
Neuroinformatics 8:90.
http://dx.doi.org/10.3389/fninf.2014.00090