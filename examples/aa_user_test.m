% Automatic analysis
% User master script based on
% github.com/rhodricusack/automaticanalysis/wiki/Manual:
% Example (aa version 4.*)
%
% Tibor Auer, MRC-CBSU
% 08-05-2013

%% INITIALISE
clear

aa_ver4_nocloud

%% DEFINE SPECIFIC PARAMETERS
% ANALYSIS RECIPE
%  Default recipe without model
aap=aarecipe('aap_parameters_defaults.xml','aap_tasklist_fmri.xml');

% Modify standard recipe module selection here if you'd like
aap.options.wheretoprocess = 'qsub'; % parallel; typical localsingle
aap.options.NIFTI4D = 1; % 4D support; typical value 0
aap.options.autoidentifyfieldmaps=1;  % required for unwarping; typical value 1
aap.options.email='All.Knowing@mrc-cbu.cam.ac.uk';

% STUDY
% Directory for analysed data
aap.acq_details.root = '/imaging/xy00/World_Universe_and_Everything';
aap.directory_conventions.analysisid = 'Nature_Paper';

aap=aas_addsubject(aap,mri_findvol(90952),[7]);
aap = aas_addsession(aap,'Loc');

% Obtain TR from the first session
h = dicominfo(mri_finddcm(90952,7));
TR = h.RepetitionTime/1000; % in seconds

aap = aas_addevent(aap,'aamod_firstlevel_model','*','*','RightFinger',[70.130 140.288 170.339 240.480 310.605 340.671]-aap.acq_details.numdummies*TR,[10 10 10 10 10 10]);
aap = aas_addevent(aap,'aamod_firstlevel_model','*','*','Rest',[30.057 100.198 160.318 200.390 270.531 330.649]-aap.acq_details.numdummies*TR,[10 10 10 10 10 10]);
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts','*','singlesession: Loc',[1 -1],'Loc_RightFinger-Rest','T');

%% DO ANALYSIS
aa_doprocessing(aap);
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
aas_garbagecollection(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid),true);
% clear all;