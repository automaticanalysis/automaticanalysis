clear;
aa_ver4

%% RECIPE:
aap = aarecipe('aap_parameters_defaults_CBSU.xml','aap_tasklist_meg.xml');

% SITE-SPECIFIC CONFIGURATION:
aap.options.wheretoprocess='qsub'; % queuing system			% typical value localsingle or qsub
aap.options.email='All.Knowing@mrc-cbu.cam.ac.uk';
aap.tasksettings.aamod_meg_denoise_ICA_2_applytrajectory.toremove = 'spat';
aap.tasksettings.aamod_meg_epochs.timewindow = [-2000 500];

%% DATA
% Directory & sub-directory for analysed data:
aap.acq_details.root = '/imaging/xy00/World_Universe_and_Everything'; 
aap.directory_conventions.analysisid = 'Nature_Paper'; 

% Add extra files
aap= aas_addinitialstream(aap,'channellabels',{'/imaging/rh01/VectorView_MAG_GRD_EEG_EOG_STI101.mat'});
aap= aas_addinitialstream(aap,'topography',{'/imaging/rh01/Methods/MEGEEGArtifactTemplateTopographies.mat'});

% Directory for raw data:
aap.directory_conventions.rawmegdatadir = '/megdata/cbu/ftd';

% Add subject (full):
aap = aas_add_meg_session(aap,'psp_bp');
aap = aas_addsubject(aap,{[12 442] []},{'psp_button_press_self_raw.fif'});
aap = aas_addsubject(aap,{[13 133] []},{'ftd_0133_bsp_raw.fif'});

% Add conditions
aap = aas_add_meg_event(aap,'aamod_meg_epochs',[12 442],'psp_bp','BP',{'STI101_down' 8192},34);
aap = aas_add_meg_event(aap,'aamod_meg_epochs',[13 133],'psp_bp','BP',{'STI101_down' 22},34);

%% RUN
aa_doprocessing(aap);
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));