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
     'CC320500' 131117; ...
     };

aa_ver5

%% DEFINE SPECIFIC PARAMETERS
%  Default recipe without model
aap=aarecipe('aap_parameters_defaults_CBSU.xml','aap_tasklist_ASL.xml');
aap = aas_configforSPM12(aap);

% Modify standard recipe module selection here if you'd like
aap.options.wheretoprocess = 'qsub'; % queuing system	% typical value localsingle or qsub
aap.options.NIFTI4D = 1;
% aap.options.email='tibor.auer@mrc-cbu.cam.ac.uk';

aap.options.autoidentifystructural_chooselast = 1;

aap=aas_renamestream(aap,'aamod_maths_00001','input','ASL_CBF');
aap.tasksettings.aamod_maths.operation='div(X,10)'; % ml/100g/min
aap=aas_renamestream(aap,'aamod_maths_00001','output','ASL_CBF','output');

%% STUDY
% Directory for analysed data
aap.acq_details.root = '/imaging/ta02/aa'; 
aap.directory_conventions.analysisid = 'ASL'; 

% Add data
aap.directory_conventions.rawdatadir = '/mridata/camcan280';
aap=aas_add_special_session(aap,'ASL');
for s = 1:size(SUBJ,1)
   ser = [...
       sscanf(basename(spm_select('FPListRec',mri_findvol(aap,SUBJ{s,2},1),'dir','.*_ep2d_tra_pasl$')),aap.directory_conventions.seriesoutputformat) ...
       sscanf(basename(spm_select('FPListRec',mri_findvol(aap,SUBJ{s,2},1),'dir','.*_Perfusion_Weighted$')),aap.directory_conventions.seriesoutputformat) ...
       sscanf(basename(spm_select('FPListRec',mri_findvol(aap,SUBJ{s,2},1),'dir','.*_relCBF$')),aap.directory_conventions.seriesoutputformat) ...
       ];
   aap = aas_addsubject(aap,SUBJ{s,1},SUBJ{s,2},'specialseries',{ser});
end

%% DO ANALYSIS
aa_doprocessing(aap);