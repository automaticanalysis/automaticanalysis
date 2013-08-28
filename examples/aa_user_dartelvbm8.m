% AA_USER_DARTELVBM8 is an AA user script illustrating the preprocessing
% for a voxel-based morphometric analysis, using the "new" segmentation
% introduced in SPM8 and DARTEL registration routines.
%
% As with all AA user scripts, the actual stages of the analysis are
% specified in an XML file - here, in aap_tasklist_dartelvbm8.xml. Although
% pre-defined recipes are distributed with AA, it is easy to modify these
% to suit your needs. You can either save a copy of the same file in your
% matlab path such that it overrides the built-in recipe, or (probably
% more easily) simply save a copy with a new name.
%
% Although the recipe file will set up the stages of analysis and default
% options, additional options can be set for both the overall analysis and
% specific modules.

clear all

%% Initialize AA and set basic recipe (aap_tasklist_dartelvbm8.xml)

aa_ver4_nocloud()
aap = aarecipe('aap_parameters_defaults.xml','aap_tasklist_dartelvbm8.xml');

  
%% Set options
    
% The name of your scanner
%aap.directory_conventions.rawdatadir='Machine_MRC-CBU_MRC35119_TrioTim';
aap.acq_details.root = '/imaging/jp01/aa4_vbm_example/subj'; % the study directory
aap.options.autoidentifyfieldmaps = 0;
aap.options.autoidentifystructural_chooselast = 1;
aap.options.aa_minver = 4;
aap.directory_conventions.remotefilesystem = 'none';         % none, s3
aap.options.wheretoprocess = 'localsingle';                  % aws, localsingle, localparallel

% Set any other options
aap.tasksettings.aamod_segment8.samp = 2;

%% Add subjects
% For many subjects, this can also be scripted to loop through a list.

aap = aas_addsubject(aap, 'CBU110656_MR1/*');
aap = aas_addsubject(aap, 'CBU050015_*/*');
aap = aas_addsubject(aap, 'CBU060593_*/*');
    
    
%% Do the processing
aap = aa_doprocessing(aap);
    

