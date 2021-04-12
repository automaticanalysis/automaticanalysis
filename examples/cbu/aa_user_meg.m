% Automatic analysis
% User master script example (aa version 5.*.*)
%
% This script demonstrates a basic MEG pipeline.
%
% For internal use at MRC CBU, Cambridge, UK - requires access to the CBU imaging
% system.
%
% v2: Johan Carlin, MRC CBU, 08-08-2018
% v1: Tibor Auer, MRC-CBSU, 01-02-2016

clear;
aa_ver5

%% RECIPE:
aap = aarecipe('aap_parameters_defaults_CBSU.xml','aap_tasklist_meg.xml');
SPM = aas_inittoolbox(aap,'spm');
SPM.load;

% SITE-SPECIFIC CONFIGURATION:
aap.options.wheretoprocess='qsub'; % queuing system			% typical value localsingle or qsub

%% PIPELINE
aap.tasksettings.aamod_meeg_denoise_ICA_2_applytrajectory.toremove = 'spat';
aap.tasksettings.aamod_meeg_epochs_meg.timewindow = [-2000 500];

%% DATA
% Directory & sub-directory for analysed data:
aap.acq_details.root = fullfile(aap.acq_details.root,'aa_demo');
aap.directory_conventions.analysisid = 'meg'; 

% Add extra files
aap= aas_addinitialstream(aap,'channellabels',...
    {'/imaging/local/software/AA/test_resources/meg/VectorView_MAG_GRD_EEG_EOG_STI101.mat'});
aap= aas_addinitialstream(aap,'topography',...
    {'/imaging/local/software/AA/test_resources/meg/MEGEEGArtifactTemplateTopographies.mat'});

% Directory for raw data:
aap.directory_conventions.rawmeegdatadir = '/megdata/cbu/ftd';
aap.directory_conventions.meegsubjectoutputformat = 'meg%02d_%04d*';
aap.directory_conventions.subject_directory_format = 3;

% Add subject (full):
aap = aas_add_meeg_session(aap,'run1');
aap = aas_addsubject(aap,'S1',{[12 442] []},'functional',{'psp_button_press_self_raw.fif'});
aap = aas_addsubject(aap,'S2',{[13 133] []},'functional',{'ftd_0133_bsp_raw.fif'});

% Add conditions
aap = aas_add_meeg_event(aap,'aamod_meeg_epochs_meg','S1','run1','BP',{'STI101_down' 8192},34);
aap = aas_add_meeg_event(aap,'aamod_meeg_epochs_meg','S2','run1','BP',{'STI101_down' 22},34);

%% RUN
aa_doprocessing(aap);
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
