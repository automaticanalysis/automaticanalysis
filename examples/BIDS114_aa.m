% This is the user master script to process BIDS multimodal NIfTI dataset ds114 (https://github.com/INCF/BIDS-examples/tree/master/ds114)

%% INITIALISE
clear

aa_ver4

%% DEFINE SPECIFIC PARAMETERS
%  Default recipe without model
aap=aarecipe('aap_parameters_defaults_CBSU.xml','BIDS114_tasklist.xml');
aap = aas_configforSPM12(aap);

% Modify standard recipe module selection here if you'd like
aap.options.wheretoprocess = 'qsub'; % queuing system			% typical value qsub | localsingle
aap.options.NIFTI4D = 1;										% typical value 0
aap.options.email='xy01@mrc-cbu.cam.ac.uk';

aap.tasksettings.aamod_dartel_norm_write.vox = 1;
aap.tasksettings.aamod_diffusion_bet.bet_f_parameter = 0.4;
aap.tasksettings.aamod_slicetiming.autodetectSO = 1;
aap.tasksettings.aamod_slicetiming.refslice = 16;
aap.tasksettings.aamod_norm_write_dartel.vox = [3 3 3];
aap.tasksettings.aamod_smooth.FWHM = 5;
aap.tasksettings.aamod_firstlevel_model(1).includemovementpars = 1; % Include/exclude Moco params in/from DM, typical value 1
aap.tasksettings.aamod_firstlevel_model(2).includemovementpars = 1; % Include/exclude Moco params in/from DM, typical value 1
aap.tasksettings.aamod_firstlevel_threshold(1).overlay.nth_slice = 9;
aap.tasksettings.aamod_firstlevel_threshold(2).overlay.nth_slice = 9;

%% STUDY
% Directory for analysed data
aap.acq_details.root = '/imaging/xy01/aa'; 
aap.directory_conventions.analysisid = 'BIDS114'; 

% Add data
aap.directory_conventions.rawdatadir = '/Data/BIDS/ds114';
aap.acq_details.numdummies = 0;
aap.options.autoidentifystructural_choosefirst = 1;
aap = aas_processBIDS(aap);

aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','singlesession:overtverbgeneration',[1],'Task','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','singlesession:overtwordrepetition',[1],'Task','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','sameforallsessions',[1],'Task','T');
 
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:covertverbgeneration',[1],'Task','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:fingerfootlips',[1 0 0],'Finger','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:fingerfootlips',[0 1 0],'Foot','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:fingerfootlips',[0 0 1],'Lips','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:linebisection',[1 0 0 -1 0],'Task:Resp-NoResp','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:linebisection',[-1 0 0 1 0],'Task:NoResp-Resp','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:linebisection',[0 0 -1 0 1],'Control:Resp-NoResp','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:linebisection',[0 0 1 0 -1],'Control:NoResp-Resp','T');

%% DO ANALYSIS
aa_doprocessing(aap);
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));