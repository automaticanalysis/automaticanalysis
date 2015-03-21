% Automatic analysis
% User master script based on
% github.com/rhodricusack/automaticanalysis/wiki/Manual:
% Example (aa version 4.*)
%
% Tibor Auer, MRC-CBSU
% 09-12-2013

%% INITIALISE
clear

aa_ver4_nocloud

%% DEFINE SPECIFIC PARAMETERS
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
aap.tasksettings.aamod_firstlevel_model.xBF.UNITS = 'secs';        	% OPTIONS: 'scans'|'secs' for onsets and durations, typical value 'secs'
aap.tasksettings.aamod_firstlevel_model.includemovementpars = 0;% Include/exclude Moco params in/from DM, typical value 1

%% STUDY
% Directory for analysed data
aap.acq_details.root = '/imaging/xy00/World_Universe_and_Everything'; 
aap.directory_conventions.analysisid = 'Nature_Paper'; 

% Add data
aap = aas_addsession(aap,'Loc');
aap = aas_addsubject(aap,90973,[7]);
aap = aas_addsubject(aap,90979,[7]);

% Add model
% Obtain TR from the first session
h = spm_dicom_headers(mri_finddcm(aap, 90973,7));
TR = h{1}.RepetitionTime/1000; % in seconds

aap = aas_addevent(aap,'aamod_firstlevel_model',mri_findvol(aap,90973),'*',...
    'REST',...                                                                                  % name
    [30.0560 100.1620 160.2800 200.3700 280.5340 340.6670]-aap.acq_details.numdummies*TR,...    % onsets
    [10.0070  10.0230  10.0230  10.0220  10.0210  10.0230]);                                    % durations
aap = aas_addevent(aap,'aamod_firstlevel_model',mri_findvol(aap,90973),'*',...
    'RIGHTFINGER',...                                                                           % name
    [70.0960 140.2360 170.3030 240.4600 320.6220 350.6890]-aap.acq_details.numdummies*TR,...    % onsets
    [10.0220  10.0220  10.0220  10.0220  10.0230  10.0230]);                                    % durations

aap = aas_addevent(aap,'aamod_firstlevel_model',mri_findvol(aap,90979),'*',...
    'REST',...                                                                                  % name
    [30.0570 100.1310 160.1990 200.2380 270.3610 330.4300]-aap.acq_details.numdummies*TR,...    % onsets
    [10.0070  10.0220  10.0060  10.0220  10.0070  10.0220]);                                    % durations
aap = aas_addevent(aap,'aamod_firstlevel_model',mri_findvol(aap,90979),'*',...
    'RIGHTFINGER',...                                                                           % name
    [70.0970 140.1860 170.2050 240.3110 310.4020 340.4510]-aap.acq_details.numdummies*TR,...    % onsets
    [10.0060  10.0070  10.0210  10.0220  10.0210  10.0220]);                                    % durations

aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts','*','singlesession:Loc',[-1 1],'Loc_RightFinger-Rest','T');

%% DO ANALYSIS
aa_doprocessing(aap);
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
aas_garbagecollection(aap,true);