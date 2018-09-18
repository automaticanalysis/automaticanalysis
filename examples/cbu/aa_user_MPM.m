% Automatic analysis (aa) - user master script
% 
% This script demonstrates multi-parameter mapping (MPM) using the CamCAN-frail data.
%
% For internal use at MRC CBU, Cambridge, UK - requires access to the CBU imaging
% system.
%
% v2: Johan Carlin, MRC CBU, 08-08-2018
% v1: Tibor Auer, MRC-CBSU, 09-12-2013

%% INITIALISE
clear
aa_ver5
 

%% DEFINE STUDY SPECIFIC PARAMETERS
aap=aarecipe('aap_tasklist_MPM.xml');

SUBJ = {...
     'S01' 160886;...
     };

% Modify standard recipe module selection here if you'd like
aap.options.wheretoprocess = 'qsub'; % queuing system	% typical value localsingle or qsub
aap.options.autoidentifyfieldmaps = 1;

aap.tasksettings.aamod_convert_specialseries.NIFTI4D = 0;

% Directory for analysed data
aap.acq_details.root = fullfile(aap.acq_details.root,'aa_demo');
aap.directory_conventions.analysisid = 'mpm'; 

% Add data
aap.directory_conventions.protocol_fieldmap= 'gre_field_mapping';
aap.directory_conventions.subject_directory_format = 3;
aap.directory_conventions.rawdatadir = '/mridata/camcan_f';
aap = aas_add_special_session(aap,'MPM');
for s = 1:size(SUBJ,1)
    % [B1, MT, PD, T1]
    ser = [...
        sscanf(basename(spm_select('FPListRec',mri_findvol(aap,SUBJ{s,2},1),'dir','.*_al_B1mapping_.*')),aap.directory_conventions.seriesoutputformat) ...
        sscanf(basename(spm_select('FPListRec',mri_findvol(aap,SUBJ{s,2},1),'dir','.*_mt_al_mtflash3d_.*'))',aap.directory_conventions.seriesoutputformat) ...
        sscanf(basename(spm_select('FPListRec',mri_findvol(aap,SUBJ{s,2},1),'dir','.*_pd_al_mtflash3d_.*'))',aap.directory_conventions.seriesoutputformat) ...
        sscanf(basename(spm_select('FPListRec',mri_findvol(aap,SUBJ{s,2},1),'dir','.*_t1_al_mtflash3d_.*'))',aap.directory_conventions.seriesoutputformat) ...
        ];
    aap = aas_addsubject(aap,SUBJ{s,1},SUBJ{s,2},'specialseries',{ser});
end

%% DO ANALYSIS
aa_doprocessing(aap);
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
