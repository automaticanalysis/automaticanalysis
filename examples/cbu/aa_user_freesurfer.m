% Automatic analysis
% User master script example (aa version 5.*.*)
%
% This script is a minimal demonstration of a freesurfer recon all pipeline. For a more
% complete example that includes using the surfaces to visualise effects, see
% aa_user_fmri_advanced.
%
% For internal use at MRC CBU, Cambridge, UK - requires access to the CBU imaging
% system.
%
% v2: Johan Carlin, MRC CBU, 08-08-2018
% v1: Tibor Auer, MRC-CBSU, 01-02-2016

%% INITIALISE
clear
aa_ver5

%% DEFINE SPECIFIC PARAMETERS
aap=aarecipe('aap_tasklist_freesurfer.xml');

% Modify standard recipe module selection here if you'd like
aap.options.wheretoprocess = 'qsub';

%% STUDY
% Directory for analysed data
aap.directory_conventions.subject_directory_format = 1;
aap.acq_details.root = fullfile(aap.acq_details.root,'aa_demo');
aap.directory_conventions.analysisid = 'freesurfer'; 

aap.directory_conventions.continueanalysis = 1;

% cbu
aap = aas_addsubject(aap,90952); 
aap = aas_addsubject(aap,90971);
% camcan
aap = aas_addsubject(aap,110220);
aap = aas_addsubject(aap,110252);

%% DO ANALYSIS
aa_doprocessing(aap);
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
