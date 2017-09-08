% Automatic analysis
% User master script based on
% github.com/rhodricusack/automaticanalysis/wiki/Manual:
% Example (aa version 5.*)
%
% Tibor Auer, MRC-CBSU
% 09-12-2013

%% INITIALISE
clear

SUBJ = {...
     'S01' 160886;...
     };
 
aa_ver5

%% DEFINE SPECIFIC PARAMETERS
%  Default recipe without model
aap=aarecipe('aap_parameters_defaults_CBSU.xml','aap_tasklist_MPM.xml');
aap = aas_configforSPM12(aap);

% Modify standard recipe module selection here if you'd like
aap.options.wheretoprocess = 'qsub'; % queuing system	% typical value localsingle or qsub
aap.options.NIFTI4D = 1;
aap.options.email='xy00@mrc-cbu.cam.ac.uk';

aap.tasksettings.aamod_convert_specialseries.numdummies = 0;
aap.tasksettings.aamod_convert_specialseries.NIFTI4D = 0;

% Directory for analysed data
aap.acq_details.root = '/imaging/xy00/aa'; 
aas_makedir(aap,aap.acq_details.root);
aap.directory_conventions.analysisid = 'MPM';
aap.directory_conventions.subject_directory_format = 3;

% Add data
aap.directory_conventions.rawdatadir = '/mridata/camcan_f';
aap = aas_add_special_session(aap,'MPM');
for s = 1:size(SUBJ,1)
    ser = [...
        sscanf(basename(spm_select('FPListRec',mri_findvol(aap,SUBJ{s,2},1),'dir','.*_gre_field_mapping_.*'))',aap.directory_conventions.seriesoutputformat) ...    
        sscanf(basename(spm_select('FPListRec',mri_findvol(aap,SUBJ{s,2},1),'dir','.*_al_B1mapping_.*')),aap.directory_conventions.seriesoutputformat) ...
        sscanf(basename(spm_select('FPListRec',mri_findvol(aap,SUBJ{s,2},1),'dir','.*_mt_al_mtflash3d_.*'))',aap.directory_conventions.seriesoutputformat) ...
        sscanf(basename(spm_select('FPListRec',mri_findvol(aap,SUBJ{s,2},1),'dir','.*_pd_al_mtflash3d_.*'))',aap.directory_conventions.seriesoutputformat) ...
        sscanf(basename(spm_select('FPListRec',mri_findvol(aap,SUBJ{s,2},1),'dir','.*_t1_al_mtflash3d_.*'))',aap.directory_conventions.seriesoutputformat) ...
        ];
    % [B0 mag, B0 phase, B1, MT mag, PD mag, T1 mag]
    aap = aas_addsubject(aap,SUBJ{s,1},SUBJ{s,2},'specialseries',{[ser(1) ser(1)+1 ser(2) ser(3) ser(4) ser(5)]});
end

%% DO ANALYSIS
aa_doprocessing(aap);
