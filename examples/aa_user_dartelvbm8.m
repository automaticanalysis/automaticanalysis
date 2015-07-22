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

aa_ver4
aap = aarecipe('aap_parameters_defaults_CBSU.xml','aap_tasklist_dartelvbm8.xml');

  
%% Set options
aap.options.autoidentifyfieldmaps = 0;
aap.options.autoidentifystructural_chooselast = 1;
aap.options.wheretoprocess = 'qsub'; % parallel; typical localsingle
aap.options.email='All.Knowing@mrc-cbu.cam.ac.uk';

% Set any other options
aap.tasksettings.aamod_segment8.samp = 2;

%% Study directory
aap.acq_details.root = '/imaging/xy00/World_Universe_and_Everything';
aap.directory_conventions.analysisid = 'Nature_Paper';

%% Add subjects
% For many subjects, this can also be scripted to loop through a list.
% cbu
aap = aas_addsubject(aap,90952); 
aap = aas_addsubject(aap,90971);
% camcan
aap = aas_addsubject(aap,110220);
aap = aas_addsubject(aap,110252);
    
%% Do the processing
aap = aa_doprocessing(aap);
aas_garbagecollection(aap,true);
    

