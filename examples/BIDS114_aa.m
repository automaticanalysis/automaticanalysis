% Automatic analysis
% User master script example (aa version 5.*.*)
%
% This is an example how to process BIDS multimodal NIfTI dataset "ds114" 
% (https://github.com/INCF/BIDS-examples/tree/master/ds114
% It uses a tasklist BIDS114_tasklist.xml (included)
%
% Tibor Auer, MRC-CBSU
% 08-02-2016

%% INITIALISE
clear

aa_ver5

%% DEFINE SPECIFIC PARAMETERS
%  Default recipe without model
aap=aarecipe('aap_parameters_defaults_CBSU.xml','BIDS114_tasklist.xml');
aap = aas_configforSPM12(aap);

% Modify standard recipe module selection here if you'd like
aap.options.wheretoprocess = 'qsub'; % queuing system			% typical value qsub | localsingle
aap.options.NIFTI4D = 1;										% typical value 0
aap.options.email='xy00@mrc-cbu.cam.ac.uk';

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
aap.acq_details.root = '/imaging/xy00/aa'; 
aap.directory_conventions.analysisid = 'BIDS114'; 

% Add data
aap.directory_conventions.rawdatadir = '/Data/BIDS/ds114';
aap.acq_details.numdummies = 0;
aap.acq_details.input.combinemultiple = 1;
aap.options.autoidentifystructural_choosefirst = 1;
aap = aas_processBIDS(aap);

aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','singlesession:overtverbgeneration_test',[1],'OVG_T:Task','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','singlesession:overtwordrepetition_test',[1],'OWR_T:Task','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','singlesession:overtverbgeneration_retest',[1],'OVG_R:Task','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','singlesession:overtwordrepetition_retest',[1],'OWR_R:Task','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','sameforallsessions',[1],'OAll:Task','T');
 
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:covertverbgeneration_test',[1],'CVG_T:Task','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:fingerfootlips_test',[1 0 0],'Loc_T:Finger','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:fingerfootlips_test',[0 1 0],'Loc_T:Foot','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:fingerfootlips_test',[0 0 1],'Loc_T:Lips','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:linebisection_test',[1 0 0 -1 0],'LB_T:Task:Resp-NoResp','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:linebisection_test',[-1 0 0 1 0],'LB_T:Task:NoResp-Resp','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:linebisection_test',[0 0 -1 0 1],'LB_T:Control:Resp-NoResp','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:linebisection_test',[0 0 1 0 -1],'LB_T:Control:NoResp-Resp','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:covertverbgeneration_retest',[1],'CVG_R:Task','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:fingerfootlips_retest',[1 0 0],'Loc_RT:Finger','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:fingerfootlips_retest',[0 1 0],'Loc_T:Foot','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:fingerfootlips_retest',[0 0 1],'Loc_T:Lips','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:linebisection_retest',[1 0 0 -1 0],'LB_RT:Task:Resp-NoResp','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:linebisection_retest',[-1 0 0 1 0],'LB_RT:Task:NoResp-Resp','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:linebisection_retest',[0 0 -1 0 1],'LB_RT:Control:Resp-NoResp','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:linebisection_retest',[0 0 1 0 -1],'LB_RT:Control:NoResp-Resp','T');

%% DO ANALYSIS
aa_doprocessing(aap);
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));