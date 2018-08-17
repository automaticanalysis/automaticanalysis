% Automatic analysis (aa) - user master script
%
% This script demonstrates an up-to-date SPM12 fMRI pipeline - realign-unwarp for joint
% motion correction and undistortion, segment8 template normalisation, robust
% coregistration through an initial rigid-body step.
%
% For internal use at MRC CBU, Cambridge, UK - requires access to the CBU imaging
% system.
%
% v2: Johan Carlin, MRC CBU, 08-08-2018
% v1: Tibor Auer, MRC-CBSU, 01-02-2016

%% INITIALISE
clear
aa_ver5

%% DEFINE SPECIFIC PARAMETERS
%  Default recipe with model
aap=aarecipe('aap_tasklist_fmri.xml');

% Modify standard recipe module selection here if you'd like
aap.options.wheretoprocess = 'qsub'; %'qsub'; % queuing system			% typical value localsingle
aap.options.autoidentifyfieldmaps = 1;
% Set slice order for slice timing correction
aap.tasksettings.aamod_slicetiming.autodetectSO = 1;
aap.tasksettings.aamod_slicetiming.refslice = 16;              	% reference slice (first acquired)
aap.tasksettings.aamod_norm_write.vox = [3 3 3];
aap.tasksettings.aamod_norm_write_meanepi.vox = [3 3 3];

aap.tasksettings.aamod_firstlevel_model.xBF.UNITS = 'secs';        	% OPTIONS: 'scans'|'secs' for onsets and durations, typical value 'secs'
aap.tasksettings.aamod_firstlevel_model.includemovementpars = 0;% Include/exclude Moco params in/from DM, typical value 1
aap.tasksettings.aamod_firstlevel_threshold.threshold.p = 0.01;

aap = aas_renamestream(aap,'aamod_secondlevel_threshold_00001','structural','normalised_structural');
aap.tasksettings.aamod_secondlevel_threshold.threshold.correction = 'none';

%% STUDY
% Directory for analysed data
aap.acq_details.root = fullfile(aap.acq_details.root,'aa_demo');
aap.directory_conventions.analysisid = 'fmri'; 

% Add data
aap.directory_conventions.subject_directory_format = 1;
aap = aas_addsession(aap,'Loc');
aap = aas_addsubject(aap,90973,'functional',[7]);
aap = aas_addsubject(aap,90979,'functional',[7]);

% Add model
% Obtain TR from the first session
h = spm_dicom_headers(mri_finddcm(aap, 90973,7));
TR = h{1}.RepetitionTime/1000; % in seconds

aap = aas_addevent(aap,'aamod_firstlevel_model',aas_getsubjname(aap,1),'*',...
    'REST',...                                                                                  % name
    [30.0560 100.1620 160.2800 200.3700 280.5340 340.6670]-aap.acq_details.numdummies*TR,...    % onsets
    [10.0070  10.0230  10.0230  10.0220  10.0210  10.0230]);                                    % durations
aap = aas_addevent(aap,'aamod_firstlevel_model',aas_getsubjname(aap,1),'*',...
    'RIGHTFINGER',...                                                                           % name
    [70.0960 140.2360 170.3030 240.4600 320.6220 350.6890]-aap.acq_details.numdummies*TR,...    % onsets
    [10.0220  10.0220  10.0220  10.0220  10.0230  10.0230]);                                    % durations

aap = aas_addevent(aap,'aamod_firstlevel_model',aas_getsubjname(aap,2),'*',...
    'REST',...                                                                                  % name
    [30.0570 100.1310 160.1990 200.2380 270.3610 330.4300]-aap.acq_details.numdummies*TR,...    % onsets
    [10.0070  10.0220  10.0060  10.0220  10.0070  10.0220]);                                    % durations
aap = aas_addevent(aap,'aamod_firstlevel_model',aas_getsubjname(aap,2),'*',...
    'RIGHTFINGER',...                                                                           % name
    [70.0970 140.1860 170.2050 240.3110 310.4020 340.4510]-aap.acq_details.numdummies*TR,...    % onsets
    [10.0060  10.0070  10.0210  10.0220  10.0210  10.0220]);                                    % durations

aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts','*','singlesession:Loc',[-1 1],'Loc_RightFinger-Rest','T');

%% DO ANALYSIS
aa_doprocessing(aap);
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
