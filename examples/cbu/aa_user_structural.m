% Automatic analysis
% User master script example (aa version 5.*.*)
%
% This script demonstrates a structural pipeline. The results it generates are 
% necessary inputs for aa_user_fmri_connect.
%
% For internal use at MRC CBU, Cambridge, UK - requires access to the CBU imaging
% system.
%
% v2: Johan Carlin, MRC CBU, 08-08-2018
% v1: Tibor Auer, MRC-CBSU, 01-02-2016

%% INITIALISE
clear
aa_ver5

SUBJ = {...
     'S01' 140905; ...
     'S02' 140910; ...
     'S03' 140913; ...
     'S04' 140928; ...
     'S05' 140931; ...
     };


%% DEFINE SPECIFIC PARAMETERS
%  Default recipe without model
aap=aarecipe('aap_tasklist_structural.xml');

% Modify standard recipe module selection here if you'd like
aap.options.wheretoprocess = 'qsub'; % queuing system	% typical value localsingle or qsub
aap.options.autoidentifyfieldmaps = 1;
aap.options.autoidentifystructural_chooselast = 1;
aap.options.autoidentifyt2=1;
aap.options.autoidentifyt2_chooselast = 1;
aap.tasksettings.aamod_segment8_multichan.writenormimg=0;
aap.tasksettings.aamod_dartel_normmni.fwhm=1;

%% STUDY
% Directory for analysed data
aap.acq_details.root = fullfile(aap.acq_details.root,'aa_demo');
aap.directory_conventions.analysisid = 'structural'; 

% Add data
aap.directory_conventions.subject_directory_format = 3;
for s = 1:size(SUBJ,1)
   aap = aas_addsubject(aap,SUBJ{s,1},SUBJ{s,2});
end

%% DO ANALYSIS
aa_doprocessing(aap);
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
