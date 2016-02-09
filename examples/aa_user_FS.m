% Automatic analysis
% User master script example (aa version 5.*.*)
%
% Tibor Auer, MRC-CBSU
% 01-02-2016

%% INITIALISE
clear

aa_ver5

%% DEFINE SPECIFIC PARAMETERS
% ANALYSIS RECIPE
%  Default recipe without model
aap=aarecipe('aap_parameters_defaults_CBSU.xml','aap_tasklist_freesurfer.xml');

% Modify standard recipe module selection here if you'd like
aap.options.wheretoprocess = 'qsub';
aap.options.autoidentifyfieldmaps=0;  % typical value 1
aap.options.email='All.Knowing@mrc-cbu.cam.ac.uk';

%% STUDY
% Directory for analysed data
aap.acq_details.root = '/imaging/xy00/World_Universe_and_Everything'; 
aap.directory_conventions.analysisid = 'Nature_Paper'; 

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