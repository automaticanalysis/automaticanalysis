
% Automatic analysis
% User master example (aa version 5.*.*)
%
% Tibor Auer, MRC-CBSU
% 01-02-2016

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
aap=aarecipe('aap_parameters_defaults_CBSU.xml','aap_tasklist_diffusion.xml');
aap = aas_renamestream(aap,'aamod_diffusion_dartel_denormDKI_00001','normalised_grey','normalised_white');
aap = aas_renamestream(aap,'aamod_diffusion_dartel_denormDKI_00001','native_grey','native_white','output');

aap = aas_configforSPM12(aap);

% Modify standard recipe module selection here if you'd like
aap.options.wheretoprocess = 'qsub'; % queuing system	% typical value localsingle or qsub
aap.options.NIFTI4D = 1;
aap.options.email='xy00@mrc-cbu.cam.ac.uk';

aap.options.autoidentifystructural_chooselast = 1;
aap.tasksettings.aamod_segment8_multichan.writenormimg=0;
aap.tasksettings.aamod_dartel_normmni.fwhm=1;

aap.tasksettings.aamod_convert_diffusion.numdummies = 0;
aap.tasksettings.aamod_diffusion_bet.bet_f_parameter = 0.2;
aap.tasksettings.aamod_diffusion_smooth.FWHM = 2.5;        
aap.tasksettings.aamod_diffusion_dartel_denormDKI.interp=4;
%% STUDY
aap=aas_addinitialstream(aap,'rois',{'/imaging/camcan/templates/craddock_ROI_841_Linda_FCpaper.nii'});

% Directory for analysed data
aap.acq_details.root = '/imaging/xy00/aa'; 
aap.directory_conventions.analysisid = 'Diffusion'; 
aap.directory_conventions.subject_directory_format = 3;

% Add data
aap = aas_add_diffusion_session(aap,'diffusion');
for s = 1:size(SUBJ,1)
   diffser = sscanf(basename(spm_select('FPListRec',mri_findvol(aap,SUBJ{s,2},1),'dir','.*_CBU_DKI_30dir_2bvals$')),aap.directory_conventions.seriesoutputformat);
   aap = aas_addsubject(aap,SUBJ{s,1},SUBJ{s,2},'diffusion',diffser);
end

%% DO ANALYSIS
aa_doprocessing(aap);