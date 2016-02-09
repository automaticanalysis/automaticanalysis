% Automatic analysis
% User master script based on
% github.com/rhodricusack/automaticanalysis/wiki/Manual:
% Example (aa version 5.*.*)
%
% Tibor Auer, MRC-CBSU
% 08-02-2016

%% INITIALISE
clear

SUBJ = {...
     'S01' 140905; ...
     'S02' 140910; ...
     'S03' 140913; ...
     'S04' 140928; ...
     'S05' 140931; ...
     };

aa_ver5

%% DEFINE SPECIFIC PARAMETERS
%  Default recipe without model
aap=aarecipe('aap_parameters_defaults_CBSU.xml','aap_tasklist_MTI.xml');
aap = aas_configforSPM12(aap);

% Modify standard recipe module selection here if you'd like
aap.options.wheretoprocess = 'qsub'; % queuing system	% typical value localsingle or qsub
aap.options.NIFTI4D = 1;
aap.options.email='xy01@mrc-cbu.cam.ac.uk';

aap.options.autoidentifystructural_chooselast = 1;
aap.tasksettings.aamod_segment8_multichan.writenormimg=0;
aap.tasksettings.aamod_dartel_norm_write.fwhm=1;

%% STUDY
% Directory for analysed data
aap.acq_details.root = '/imaging/xy01/aa'; 
aap.directory_conventions.analysisid = 'MTI'; 

% Add data
aap.acq_details.numdummies = 1;
aap=aas_add_special_session(aap,'MTI_MT');
aap=aas_add_special_session(aap,'MTI_baseline');
for s = 1:size(SUBJ,1)
   MTser = sscanf(basename(spm_select('FPListRec',mri_findvol(aap,SUBJ{s,2},1),'dir','.*_CBU_MTR_TR50_MT$')),aap.directory_conventions.seriesoutputformat);
   aap = aas_addsubject(aap,SUBJ{s,1},SUBJ{s,2},'specialseries',[MTser MTser+1]);
end

%% DO ANALYSIS
aa_doprocessing(aap);
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));