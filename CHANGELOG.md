## 4.2 (initial) ##

Described in ["the aa manuscript"](http://dx.doi.org/10.3389/fninf.2014.00090)

## 4.3.0 ##

#### Copyrights ####
  - Sorting _external_/_extrafunctions_
  - Acknowledging externals in [README.json](https://github.com/rhodricusack/automaticanalysis/blob/devel-share/external/README.json)

#### New general features ####
  - dynamic modification of streams (`aas_renamestream`)
  - selected_sessions now works also for branching

  - from NIfTI
    - 3D NIfTI inputs are fully supported
    - NIfTI fieldmaps are supported (`aamod_fieldmapfromnifti`)
    - [BIDS](bids.neuroimaging.io) datasets are supported (partial)
    - `aamod_fsl_reorienttoMNI`
	
#### New features for fMRI ####
  - New fMRI example to demonstrate some (new) features: ([aa\_user\_fmri\_advanced](https://github.com/rhodricusack/automaticanalysis/blob/devel-share/examples/aa_user_fmri_advanced.m))
    - Specifying structural series
    - Explaining Motion FingerPrint instead of simple mocoparameters in design
    - Automatic slicetiming with exact timing from DICOM header (`autodetectSO`)
    - DARTEL normalisation
    - Activation maps projected to surface using FreeSurfer (`aamod_firstlevel_threshold_register2FS`)
    - Second-level GIFT (`aamod_secondlevel_GIFT`)

  - automatic temporal modulation
  - `aamod_firstlevel_model`: option to save residual 
  - Despiking with BrainWavelet (`aamod_waveletdespike`)
  - `aamod_maths`
  - `aamod_temporalfilter`

#### New features for MEG ####
  - Diagnostics (courtesy to [Rik Henson](https://www.mrc-cbu.cam.ac.uk/people/rik.henson))
  - Update examples: [aa\_user\_meg.m](https://github.com/rhodricusack/automaticanalysis/blob/devel-share/examples/aa_user_meg.m), [aa\_user_meg\_connect](https://github.com/rhodricusack/automaticanalysis/blob/devel-share/examples/aa_user_meg_connect.m)
  - Maxfilter (`aamod_meg_maxfilt`) allows getting HPI from a session
  - Maxfilter (`aamod_meg_maxfilt`) accepts custom calibration files
  - Better delegeation of denoising tasks:
    - `aamod_meg_denoise_ICA_1`: runs ICA
    - `aamod_meg_denoise_ICA_1`: thresholds and removes
  - Epoching: `aas_add_meg_event` and `aamod_meg_epochs`
  - Averaging: `aamod_meg_average`

#### Improving robustness ####
  - Code development
    - Simplifying code
      - Setting modality --> `aas_getfile_bystream` automatically detects sessions based on modality
      - Avoid `eval`
      - Avoid `spm_jobman`
    - Removing 'orphan' functions
    - Update: spm_mods based on SPM12 r6470
				
  - Running
    - Compatibility with MATLAB pre-r2012b (local execution only)
    - `qsub` cleans jobs from previous execution
    - Iterative file retrieval (`aap.options.maximumretry`) 
	- Specifying streams as diagnostic (`aas_garbagecollection` will not touch them)
    - aaworker
      - More efficient folder handling: 
        - aaworker folders are created in ~/.aa (hidden) folder and 
        - aaworker folders are cleaned up regularly after specified time (`aap.options.aaworkercleanup`)
      - aaworker structure is passed during `qsub`

  - Summaries
    - Add default message for "Motion correction summary"
    - Rename "First level contrasts" to "First level thresholded maps"

#### Convenience/Ease of access ####
  - Code development
    - Example module template with session domain ([aamod_session](https://github.com/rhodricusack/automaticanalysis/blob/devel-share/examples/aamod_session.m))

  - Setting up
    - `aaClass` provides links to ["the aa manuscript"](http://dx.doi.org/10.3389/fninf.2014.00090), to ["the aa website"](http://automaticanalysis.org) and to examples
    - Allow specification of subject in the whole UMS in the same way: `aas_addcontrast`, `aas_addcovariate`, `aas_addevent`, `aamod_dartel_createtemplate`
    - `aap.options.checktasksettingconsistency`: Check whether settings have changed since the last execution and re-run the task accordingly regardless of the doneflag (experimental!)

  - Running
    - `qsub` error provides links to code
    - GUI for `qsub`
    - `aap.options.verbosity`:
      <br>0 - Does not show any message (still crashes in case of error)
      <br>1 - Shows errors only
      <br>2 - Shows every message (default)

  - Results
    - Vertical workflow chart
    - More consistent (and nicer :sparkles:) graph of contrasts 
    - `aamod_firstlevel_threshold` creates overlays along all three axes