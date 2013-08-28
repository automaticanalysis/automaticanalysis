% Automatic analysis
% User master script
% Rhodri Cusack MRC CBU Cambridge Jan 2006
% Danny Mitchell MRC CBU Cambridge Jul 2008

%% add paths
addpath /imaging/local/spm/spm5/
run /imaging/dm01/dev_trees/aa_svn/aa_ver3_devel
addpath /imaging/dm01/dev_trees/spm5_cbu_updates
addpath /imaging/dm01/MEG/aaMEG
addpath /imaging/dm01/MEG/aaMEG/FiffAccessToolbox

%% create aap structure from parameters and tasklist files
aap=aarecipe('aap_parameters_defaults.xml','aap_tasklist_meg_minimal.xml');
aap.options.aa_minver=1.0; % will only work on aa version 1.0 or above

%% The input directory for raw data 
aap.directory_conventions.rawmegdatadir = '/megdata/cbu/djm_cda'; 

%% The output directory for analysis 
aap.acq_details.root = '/imaging/dm01/MEG/17subs_aaTestParallel7/'; 

%% The subjects
aap=aas_addsubject(aap,{'meg07_0002',  '*CBU060321/*'});
aap=aas_addsubject(aap,{'meg07_0001djm','missing'});
aap=aas_addsubject(aap,{'meg07_0053', '*CBU080358/*'});
aap=aas_addsubject(aap,{'meg07_0054', '*CBU060074/*'});
aap=aas_addsubject(aap,{'meg07_0080', '*CBU060046/*'});
aap=aas_addsubject(aap,{'meg07_0081', '*CBU071094/*'});
aap=aas_addsubject(aap,{'meg07_0084', '*CBU070329/*'});
aap=aas_addsubject(aap,{'meg07_0085', '*CBU071125/*'});
aap=aas_addsubject(aap,{'meg07_0091', '*CBU071126/*'});
aap=aas_addsubject(aap,{'meg07_0128', '*CBU080012/*'});
aap=aas_addsubject(aap,{'meg08_0056', '*CBU080084/*'});
aap=aas_addsubject(aap,{'meg08_0059', '/imaging/dm01/MyStructurals/rdjm_MeanOf12.nii'});
aap=aas_addsubject(aap,{'meg08_0071', '*CBU080149/*'});
aap=aas_addsubject(aap,{'meg08_0074', '*CBU070809/*'});
 aap=aas_addsubject(aap,{'meg08_0075', 'missing'});
aap=aas_addsubject(aap,{'meg08_0076', '*CBU080006/*'});
aap=aas_addsubject(aap,{'meg08_0082', 'missing'});

%% The blocks
aap=aas_addsession(aap,'vstm');
aap=aas_addsession(aap,'esta');

%% The tasks
%  aap=aas_addtask(aap,'aamod_emeg_processmri',[],[]);
%  aap=aas_addtask(aap,'aamod_emeg_maxfilter',[],[]); % output prefix: N?ST?t4 (see below)
%  aap=aas_addtask(aap,'aamod_emeg_maxfilter',[],[]); % output prefix: N?ST?t2 (see below)
%  aap=aas_addtask(aap,'aamod_emeg_importfif',{'aamod_emeg_maxfilter','aamod_emeg_maxfilter_02'}); % output prefix: N?ST?t#
%  aap=aas_addtask(aap,'aamod_emeg_filter'); % output prefix: fN?ST?t# before epoching to avoid edge effects; changed defaults to [0.1 40] bandpass
%  aap=aas_addtask(aap,'aamod_emeg_ica'); % output prefix: piae_fN?ST?t# (before concatenation? - components my be different across blocks e.g. due to movement)
%  aap=aas_addtask(aap,'aamod_emeg_epoch'); % output prefix: e_fN?ST?t#
%  aap=aas_addtask(aap,'aamod_emeg_concatenate'); % output prefix: cpiae_fN?ST?t#
%  aap=aas_addtask(aap,'aamod_emeg_artefact'); % output prefix: ae_fN?ST?t# (before ICA? - Jason's recommendation)
%  
%  aap=aas_addtask(aap,'aamod_emeg_contrasts'); % output prefix: mcpiae_fN?ST?t#
%  aap=aas_addtask(aap,'aamod_emeg_plottopographies');
%  aap=aas_addtask(aap,'aamod_emeg_sensor2vol',[],'aamod_emeg_contrasts');
%  
%  aap=aas_addtask(aap,'aamod_emeg_splitsensors',[],'aamod_emeg_artefact'); % output prefix: scpiae_fN?ST?t#
% 
% % % source analyses
%  aap=aas_addtask(aap,'aamod_emeg_forwardmodel',[],{'aamod_emeg_processmri','aamod_emeg_splitsensors'});
%  aap=aas_addtask(aap,'aamod_emeg_inversion');
%  aap=aas_addtask(aap,'aamod_emeg_sourcecontrasts');
%  aap=aas_addtask(aap,'aamod_emeg_source2vol');
%  aap=aas_addtask(aap,'aamod_emeg_plotsources');
% 
% % time-frequency analysis  
%  aap=aas_addtask(aap,'aamod_emeg_timefrequency',[],'aamod_emeg_splitsensors'); % time-frequency analysis; mt#scpiae_fN?ST?t#
% aap=aas_addtask(aap,'aamod_emeg_tf2vol')

%  % plot
%  aap=aas_addtask(aap,'aamod_emeg_ploterfs',[],'aamod_emeg_forwardmodel'); % do this after forward model so we have coregistration for ROIs
%  
% % % group analyses
  aap=aas_addtask(aap,'aamod_emeg_groupsensoranalysis',[],'aamod_emeg_sensor2vol');
 % aap=aas_addtask(aap,'aamod_emeg_groupsourceanalysis',[],'aamod_emeg_source2vol');

%% SET ANY OTHER PARAMETERS YOU WOULD LIKE TO BE DIFFERENT FROM THE DEFAULTS
aap.options.userinterface=0;
%aap.tasksettings.aamod_emeg_maxfilter(1).MaxMoveMethod=4; % to 'default', with origin from sphere fit to head shape (for sensor level analysis)
%aap.tasksettings.aamod_emeg_maxfilter(2).MaxMoveMethod=2; % to 1st block of type per subject (perhaps best for source localisation?) 
% %% 1st round of artefact rejection ignoring HEOG prior to ica 
% aap.tasksettings.aamod_emeg_artefact(1).thresholdHEOG=Inf;
% aap.tasksettings.aamod_emeg_artefact(1).thresholdHEOGmean=Inf;
% %% 2nd round of artefact rejection using HEOG after ica 
% aap.tasksettings.aamod_emeg_artefact(2).InputFilter='^pi.*BLOCK.*\.mat$';
% aap.tasksettings.aamod_emeg_artefact(2).thresholdHEOG='auto';
% aap.tasksettings.aamod_emeg_artefact(2).thresholdHEOGmean='auto';
% aap.tasksettings.aamod_emeg_ica.DataFraction=0.5;
%aap.tasksettings.aamod_emeg_forwardmodel.ForwardModel='meg_os';

% aap.tasksettings.aamod_emeg_maxfilter.Overwrite=1
% aap.tasksettings.aamod_emeg_importfif.Overwrite=1;
% aap.tasksettings.aamod_emeg_filter.Overwrite=1
% aap.tasksettings.aamod_emeg_epoch.Overwrite=1;
% aap.tasksettings.aamod_emeg_artefact.Overwrite=1;
% aap.tasksettings.aamod_emeg_ica.Overwrite=1;
% aap.tasksettings.aamod_emeg_sensor2vol.Overwrite=1;
% aap.tasksettings.aamod_emeg_groupsensoranalysis.Overwrite=1;
%aap.tasksettings.aamod_emeg_inversion.Overwrite=0;

%% do it
aa_doprocessing(aap);
%aa_doprocessing_parallel(aap);
%aa_doprocessing_parallel(aap,'continue');

%% collate any figures
%aas_emeg_report(aap);