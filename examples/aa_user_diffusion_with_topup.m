
% Automatic analysis
% User master example (aa version 5.*.*)
%
% Tibor Auer, MRC-CBSU
% 16-08-2017

%% INITIALISE
clear

SUBJ = {...
     'S01' 160018; ...
     };

aa_ver5

%% DEFINE SPECIFIC PARAMETERS
%  Default recipe without model
aap=aarecipe('aap_parameters_defaults_CBSU.xml','aap_tasklist_diffusion_with_topup.xml');
aap = aas_renamestream(aap,'aamod_diffusion_dartel_denormDKI_00001','normalised_grey','normalised_white');
aap = aas_renamestream(aap,'aamod_diffusion_dartel_denormDKI_00001','native_grey','native_white','output');

aap = aas_configforSPM12(aap);

% Modify standard recipe module selection here if you'd like
aap.options.wheretoprocess = 'qsub'; % queuing system	% typical value localsingle or qsub
aap.options.NIFTI4D = 1;
aap.options.email='xy00@mrc-cbu.cam.ac.uk';

aap.options.autoidentifystructural_choosefirst = 1;
aap.tasksettings.aamod_segment8_multichan.writenormimg=0;
aap.tasksettings.aamod_dartel_norm_write.fwhm=1;

aap.tasksettings.aamod_convert_diffusion_phaseencode_direction.numdummies = 0;
aap.tasksettings.aamod_diffusion_extractnodif.tol_b0 = 10;
aap.tasksettings.aamod_diffusion_bet.bet_f_parameter = 0.2;
aap.tasksettings.aamod_diffusion_smooth.FWHM = 2.5;        
aap.tasksettings.aamod_diffusion_dartel_denormDKI.interp=4;
%% STUDY
aap=aas_addinitialstream(aap,'rois',{'/imaging/camcan/templates/craddock_ROI_841_Linda_FCpaper.nii'});

% Directory for analysed data
aap.acq_details.root = '/imaging/xy00/aa/topup'; 
aap.directory_conventions.analysisid = 'Diffusion'; 
aap.directory_conventions.subject_directory_format = 3;

% Add data
aap = aas_add_diffusion_session(aap,'diffusion');
for s = 1:size(SUBJ,1)
    diffser = cellfun(@(x) sscanf(x,aap.directory_conventions.seriesoutputformat),...
        cellstr(basename(spm_select('FPListRec',mri_findvol(aap,SUBJ{s,2},1),'dir','.*_diff_mbep2d_.*'))))';
    aap = aas_addsubject(aap,SUBJ{s,1},SUBJ{s,2},'diffusion',{diffser});
end

%% DO ANALYSIS
aa_doprocessing(aap);
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
