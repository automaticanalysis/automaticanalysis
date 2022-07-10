% Automatic analysis (aa) - user master script
%
% This script demonstrate an Arterial Spin Labelling pipeline using CamCAN data.
% 
% v2: Johan Carlin, MRC CBU, 2018-08-08
% v1: Tibor Auer, MRC-CBSU, 08-02-2016

%% INITIALISE
clear
aa_ver5

SUBJ = {...
     'CC320500' 131117; ...
     };

%% LOAD TASKLIST
aap=aarecipe('ASL.xml');

% this example uses SPM tools in the user script, so we have to ensure SPM is
% on the path
spmhit = which('spm_spm');
spmdir = aas_gettoolboxdir(aap,'spm');
if any(spmhit)
    assert(strcmp(fileparts(spmhit), spmdir), ...
        'spm on path differs from that in aap.directory_conventions');
else
    fprintf('adding spmdir to path: %s\n', spmdir);
    SPMtool = spmClass(spmdir);
    SPMtool.load;
end

% Modify standard recipe module selection here if you'd like
aap.options.wheretoprocess = 'qsub'; % queuing system	% typical value localsingle or qsub

aap.options.autoidentifystructural_chooselast = 1;
aap.tasksettings.aamod_ASL.scaleCBF = 0.1;
aap.tasksettings.aamod_ASL_coreg_extended_2.eoptions.cost_fun = 'ecc';

%% STUDY
% Directory for analysed data
aap.acq_details.root = fullfile(aap.acq_details.root,'aa_demo');
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
% prepare report (this can be slow...)
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
