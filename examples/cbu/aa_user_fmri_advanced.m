% Automatic analysis
% User master script example (aa version 5.*.*)
%
% This script demonstrates an advanced fMRI pipeline (compare with aa_user_fmri):
%	- Specifying structural series
% 	- Motion FingerPrint instead of simple mocoparameters
% 	- DARTEL normalisation
% 	- Activation maps projected to surface (FreeSurfer)
% 	- Second-level GIFT
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
aap=aarecipe('aap_tasklist_fmri_advanced.xml');

% Modify standard recipe module selection here if you'd like
aap.options.wheretoprocess = 'qsub'; % queuing system			% OPTIONS: 'localsingle'|'qsub' for aa engine, typical value 'qsub'
aap.options.autoidentifyfieldmaps = 1;
% Set slice order for slice timing correction
aap.tasksettings.aamod_realignunwarp.mfp.run = 1;               % Motion FingerPrint, typical value 0
aap.tasksettings.aamod_slicetiming.autodetectSO = 1;           	% auto
aap.tasksettings.aamod_slicetiming.refslice = 16;              	% reference slice (first acquired)
aap.tasksettings.aamod_smooth.FWHM = 5; 						% smoothing kernel size, typical value 10
aap.tasksettings.aamod_firstlevel_model.xBF.name = 'hrf (with time and dispersion derivatives)';
aap.tasksettings.aamod_firstlevel_model.xBF.UNITS = 'secs';    	% OPTIONS: 'scans'|'secs' for onsets and durations, typical value 'secs'
aap.tasksettings.aamod_firstlevel_model.includemovementpars = 1;% Include/exclude Moco params in/from DM, typical value 1
aap.tasksettings.aamod_firstlevel_threshold.threshold.p = 0.01;

aap = aas_renamestream(aap,'aamod_secondlevel_GIFT_00001','epi','aamod_smooth_00001.epi');

%% STUDY
% Directory for analysed data
aap.acq_details.root = fullfile(aap.acq_details.root,'aa_demo_fmri_advanced');
aap.directory_conventions.analysisid = 'test_fmri_advanced'; 

% Add data
aap.directory_conventions.subject_directory_format = 3;
aap = aas_addsession(aap,'Loc');
aap = aas_addsubject(aap,'S1',90973,'structural',2,'functional',7);
aap = aas_addsubject(aap,'S2',90979,'structural',2,'functional',7);

% Add model
% Obtain TR from the first session
h = dicominfo(mri_finddcm(aap, 90973,7));
TR = h.RepetitionTime/1000; % in seconds
% The "hard"(-coded) way
aap = aas_addevent(aap,'aamod_firstlevel_model','S1','*',...
    'REST',...                                                                                  % name
    [30.0560 100.1620 160.2800 200.3700 280.5340 340.6670]-aap.acq_details.numdummies*TR,...    % onsets
    [10.0070  10.0230  10.0230  10.0220  10.0210  10.0230]);                                    % durations
aap = aas_addevent(aap,'aamod_firstlevel_model','S1','*',...
    'RIGHTFINGER',...                                                                           % name
    [70.0960 140.2360 170.3030 240.4600 320.6220 350.6890]-aap.acq_details.numdummies*TR,...    % onsets
    [10.0220  10.0220  10.0220  10.0220  10.0230  10.0230]);                                    % durations

aap = aas_addevent(aap,'aamod_firstlevel_model','S2','*',...
    'REST',...                                                                                  % name
    [30.0570 100.1310 160.1990 200.2380 270.3610 330.4300]-aap.acq_details.numdummies*TR,...    % onsets
    [10.0070  10.0220  10.0060  10.0220  10.0070  10.0220]);                                    % durations
aap = aas_addevent(aap,'aamod_firstlevel_model','S2','*',...
    'RIGHTFINGER',...                                                                           % name
    [70.0970 140.1860 170.2050 240.3110 310.4020 340.4510]-aap.acq_details.numdummies*TR,...    % onsets
    [10.0060  10.0070  10.0210  10.0220  10.0210  10.0220]);   

aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts','*','singlesession:Loc',[-1 1],'Loc_RightFinger-Rest','T');

%% DO ANALYSIS
aa_doprocessing(aap);
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
