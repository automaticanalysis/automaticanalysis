% Automatic analysis
% User master script based on
% github.com/rhodricusack/automaticanalysis/wiki/Manual:
% Example (aa version 4.*)
%
% Tibor Auer, MRC-CBSU
% 20-03-2013

%% INITIALISE
clear

aa_ver4_nocloud

%% DEFINE SPECIFIC PARAMETERS
% ANALYSIS RECIPE
%  Default recipe without model
aap=aarecipe('aap_parameters_defaults_CBSU.xml','aap_tasklist_fmri.xml');

% Modify standard recipe module selection here if you'd like
aap.options.wheretoprocess = 'qsub'; % queuing system			% typical value localsingle
aap.options.autoidentifyfieldmaps=1;  							% typical value 1
aap.options.NIFTI4D = 1;										% typical value 0
aap.options.email='All.Knowing@mrc-cbu.cam.ac.uk';
% Set slice order for slice timing correction
aap.tasksettings.aamod_slicetiming.sliceorder=[32:-1:1];       	% descending
aap.tasksettings.aamod_slicetiming.refslice = 16;              	% reference slice (first acquired)
aap.tasksettings.aamod_firstlevel_model.UNITS = 'secs';        	% OPTIONS: 'scans'|'secs' for onsets and durations, typical value 'secs'
aap.tasksettings.aamod_firstlevel_model.includemovementpars = 0;% Include/exclude Moco params in/from DM, typical value 1

% STUDY
% Directory for analysed data
aap.acq_details.root = '/imaging/xy00/World_Universe_and_Everything'; 
aap.directory_conventions.analysisid = 'Nature_Paper'; 

aap=aas_addsubject(aap,90952,[7]);
aap = aas_addsession(aap,'Loc');
aap = aas_addsubject(aap,90971,[7]);
aap = aas_addsession(aap,'Loc');

% Obtain TR from the first session
h = dicominfo(mri_finddcm(aap, 90952,7));
TR = h.RepetitionTime/1000; % in seconds

aap = aas_addevent(aap,'aamod_firstlevel_model','*','*','RightFinger',[70.130 140.288 170.339 240.480 310.605 340.671]-aap.acq_details.numdummies*TR,[10 10 10 10 10 10]);
aap = aas_addevent(aap,'aamod_firstlevel_model','*','*','Rest',[30.057 100.198 160.318 200.390 270.531 330.649]-aap.acq_details.numdummies*TR,[10 10 10 10 10 10]);
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts','*','singlesession:Loc',[1 -1],'Loc_RightFinger-Rest','T');

%% DO ANALYSIS
aa_doprocessing(aap);
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
aas_garbagecollection(aap,true);
% clear all;