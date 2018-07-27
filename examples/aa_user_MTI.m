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
% aap.options.email='xy01@mrc-cbu.cam.ac.uk';

aap.options.autoidentifystructural_chooselast = 1;
aap.tasksettings.aamod_segment8_multichan.writenormimg=0;
aap.tasksettings.aamod_dartel_norm_write.fwhm=1;

aap.tasksettings.aamod_convert_specialseries.numdummies = 0;
aap.tasksettings.aamod_convert_specialseries.NIFTI4D = 0;

%% STUDY
% Directory for analysed data - NB, you can place a symbolic link in ~/aa if you
% don't want data in your home.
% (we don't support PC, but here you'd do getenv('USERPROFILE'))
aadir = fullfile(getenv('HOME'),'aa');
if ~exist(aadir,'dir')
    success = mkdir(aadir);
    assert(success, 'creating data directory failed');
end
aap.acq_details.root = aadir;
aap.directory_conventions.analysisid = 'MTI'; 


% Add data
aap.acq_details.numdummies = 1;
aap=aas_add_special_session(aap,'MTI');
for s = 1:size(SUBJ,1)
   MTser = sscanf(basename(spm_select('FPListRec',mri_findvol(aap,SUBJ{s,2},1),'dir','.*_CBU_MTR_TR50_MT$')),aap.directory_conventions.seriesoutputformat);
   aap = aas_addsubject(aap,SUBJ{s,1},SUBJ{s,2},'specialseries',{[MTser MTser+1]});
end

aap = aas_addinitialstream(aap,'rois', ...
    {'/imaging/camcan/templates/Juelich-maxprob-thr25-2mm.nii'});

%% DO ANALYSIS
aa_doprocessing(aap);
