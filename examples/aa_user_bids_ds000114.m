% Automatic analysis
% User master script example (aa version 5.*.*)
%
% This is an example how to process a BIDS multimodal NIfTI dataset. We assume that you
% have already downloaded the dataset to your directory.
%
% Tibor Auer, MRC-CBSU
% 08-02-2016

%% INITIALISE
clear
aa_ver5

%% LOAD TASKLIST
aap = aarecipe('aap_tasklist_bids_ds000114.xml');

% Modify standard recipe module selection here if you'd like
aap.options.wheretoprocess = 'qsub'; % queuing system			% typical value qsub | localsingle

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
aap.acq_details.root = fullfile(aap.acq_details.root,'aa_demo');
aap.directory_conventions.analysisid = 'bids_ds114'; 

% Add data
% download the demo dataset (if necessary)
aap = aa_downloaddemo(aap,'https://files.osf.io/v1/resources/9q7dv/providers/osfstorage/57e549f9b83f6901d457d162');
% The bids parser only supports a single rawdatadir. Pick the one that has bids in it.
sources = strsplit(aap.directory_conventions.rawdatadir,':')';
bidsind = cell_index(sources,'ds114_test2');
assert(~any(bidsind==0), ...
    'did not find "ds114_test2" BIDS directory in aap.directory_conventions.rawdatadir')
assert(numel(bidsind)==1,...
    'multiple reference to "ds114_test2" BIDS directory in aap.directory_conventions.rawdatadir')
aap.directory_conventions.rawdatadir = sources{bidsind};

aap.acq_details.numdummies = 0;

aap.acq_details.input.combinemultiple = 1;
aap.options.autoidentifystructural_choosefirst = 1;
aap = aas_processBIDS(aap,[],[],{'sub-01','sub-02'});

aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','sameforallsessions',[1 0 0],'All:Finger','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','sameforallsessions',[0 1 0],'All:Foot','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','sameforallsessions',[0 0 1],'All:Lips','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','sessions:+1xfinger_foot_lips_test|-1xfinger_foot_lips_retest',[1 0 0],'T-RT:Finger','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','sessions:+1xfinger_foot_lips_test|-1xfinger_foot_lips_retest',[0 1 0],'T-RT:Foot','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','sessions:+1xfinger_foot_lips_test|-1xfinger_foot_lips_retest',[0 0 1],'T-RT:Lips','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','singlesession:finger_foot_lips_test',[1 0 0],'Loc_T:Finger','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','singlesession:finger_foot_lips_test',[0 1 0],'Loc_T:Foot','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','singlesession:finger_foot_lips_test',[0 0 1],'Loc_T:Lips','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','singlesession:finger_foot_lips_retest',[1 0 0],'Loc_RT:Finger','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','singlesession:finger_foot_lips_retest',[0 1 0],'Loc_T:Foot','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','singlesession:finger_foot_lips_retest',[0 0 1],'Loc_T:Lips','T');

aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','sameforallsessions',[1 0 0 -1 0],'All:Task:Resp-NoResp','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','sameforallsessions',[0 0 -1 0 1],'All:Control:Resp-NoResp','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','sessions:+1xline_bisection_test|-1xline_bisection_retest',[1 0 0 -1 0],'T-RT:Task:Resp-NoResp','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','sessions:+1xline_bisection_test|-1xline_bisection_retest',[0 0 -1 0 1],'T-RT:Control:Resp-NoResp','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:line_bisection_test',[1 0 0 -1 0],'LB_T:Task:Resp-NoResp','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:line_bisection_test',[0 0 -1 0 1],'LB_T:Control:Resp-NoResp','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:line_bisection_retest',[1 0 0 -1 0],'LB_RT:Task:Resp-NoResp','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00002','*','singlesession:line_bisection_retest',[0 0 -1 0 1],'LB_RT:Control:Resp-NoResp','T');

%% DO ANALYSIS
aa_doprocessing(aap);
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
