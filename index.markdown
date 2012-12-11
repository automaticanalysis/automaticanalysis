= What is aa, and why use it? =
Automatic analysis (aa) is a pipeline system for neuroimaging, written in Matlab. It supports SPM 5/8 and some functions from FSL.

*'''Automatic'''. Virtually automatic analysis: full group-level statistics using a [[Available-Tasklists|typical fMRI recipe]] with minimal coding
*'''Flexible control'''. Users new to neuroimaging are directed to what is essential with other settings taking sensible defaults. More experienced users can easily change a range of settings, and advanced users can modify components of the system to change any behaviour.
*'''Restartable'''. If AA stops for any reason, restarting it will make it begin at the stage where it left off.
*'''Parallel processing'''. Where multiple machines are available, AA jobs can seamlessly be distributed across them.
*'''Well behaved'''. Field maps and structurals are automatically detected. Structurals can be automatically copied to a chosen location.
*'''Minimal overhead'''. To run an analysis, you only need one short script, and access to the aa library, which is typically installed in a central location at each site as a read-only directory structure. Individual settings can then be changed, new files added or files overridden. This means that the number of scripts for each new fMRI study is kept minimal, making it easier to track down what happened at a later date.
*'''Easy to maintain''' The code components are stored in a github repository, such that changes and new components can be easily be released and updated.
*'''Record keeping'''. Unlike SPM used from the GUI, the system records all parameters used, and allows easy recreation of a dataset from the raw data at a later date.
*'''Modular, with simple interface'''. Matlab programmers can easily write new [[Manual |modules]] and incorporate them into the processing stream.
*'''Support for EEG & MEG as well as FMRI'''. Exploit the power of SPM's unique source localisation and multiple comparisons correction.

= Recent updates =
'''13/3/2012 new feature'''
* You may now set up more complex pipelines by qualifying inputs to a module in the xml (or programatically) like this: 
```xml
<inputstreams>
<stream>aamod_realign00001.epi</stream>
</inputstreams>
```
Within the module, you may refer to the stream either by its qualified name (useful if there is more than one - say epis from to different stages) or by its abbreviated name (just epi) if there is only one.

'''Jan 2012 bugfix'''
* There was a bug in aamod_firstlevel_model that will have disrupted your events if you used durations that different from onset to onset, and you didn't pre-sort your onsets (i.e., they weren't in ascending temporal order)
'''20/1/2012 features'''
* added fieldmap capability with unwarping
* extended (more robust) coregistration of EPIs & structural to MNI templates
* robust brain extraction with fsl BET
* capability for multiple structurals (e.g. MP2RAGE; T1&T2)
* automatic detection of volume TR for Siemens 3D sequences
* new diagnostic images (e.g. fsl BET, normalisation and extended coregistration)
'''10/12/2011 new feature'''
Added facility to queue jobs through [http://research.cs.wisc.edu/condor/ Condor]. Put into your user script:
```matlab
aap.options.wheretoprocess='condor';
```
'''4/12/2011 new feature'''
Raw dicom data for each subject can now be arranged in any kind of arbitrary directory structure, such as all files in one folder (like Robart's data), one file per MR acquisition (like CBU data), or in fact any other structure whatsoever. The "aamod_autoidentifyseries" module scans every dicom file provided for its series and acquisition numbers, and organizes the input using these.
'''20/3/2011 new documentation'''
[[Scripting|Style Guide, How to convert version 3 modules to version 4]]
'''20/3/2011 new documentation'''
[[Stream-Reference|new page on streams]] 
'''19/3/2011 new feature'''
Ability to [[Branched-Tasklists|select particular sessions to be run on a branch of a tasklist]]
'''17/3/2011 new feature'''
New [[Modelling-and-Contrasts|first level model and contrast]] scripts with little or no scripting
The "AA" (Automatic Analysis) software package automates the analysis of neuroimaging data using SPM. It is written in Matlab and designed to be compatible with recent versions of SPM (SPM 5 and SPM 8) and NIFTI data.

It has been developed by [http://www.cusacklab.org cusacklab.org] with contributions from others at the MRC CBU in Cambridge, UK, the Brain and Mind Institute at Western University, and the Donders Centre for Cognitive Neuroimaging

= What's new in aa version 4? =
'''More flexible pipelines'''
* Simpler to specify and so fewer opportunities to make mistakes
* Easy reordering of modules
* [[Branched-Tasklists|Arbitrarily branched pipelines]]

'''Better organisation and more efficient calculation'''
* Structured record of inputs and outputs from each module
* Results are structured, making them quicker to navigate or tidy up
* If you change some earlier stage, only the later stages that can be affected are redone
* There are more opportunities for parallel execution of modules

'''New modules'''
* Easy specification of fMRI models and contrasts with ''' minimal or no Matlab coding '''. It's as easy as adding lines like this to your user script:
```matlab
aap=aas_addevent(aap,'aamod_firstlevel_model','*','*','VisualStimulus',ons,dur);
```
* Support for multi-echo EPI data
* VBM with Dartel (from Jonathan Peelle)
* MVPaa (Multi Voxel Pattern [automatic] analysis; from Alejandro Vicente Grabovetsky)

'''More flexible computing configurations'''
* The system knows what specific data any given module needs, and can just provide it with this. This allows better use of local scratch storage, rather than relying on remote (NFS or Amazon S3) file storage.
* In-built integration with cloud computing

Take a look at other [[New-Features-In-Version-4]]

----

Note: This wiki hosts the documentation for the Automatic Analysis (AA) program. Continue onward to the AA documentation by using our navigation sidebar.