% Automatic analysis
% User master script example (aa version 5.*.*)
%
% This script demonstrates a basic EEG pipeline.

clear;
aa_ver5

SUBJS = [10002 10003 10004 10005];

%% RECIPE
aap = aarecipe('aap_tasklist_meeg.xml');
SPM = spmClass(aap.directory_conventions.toolboxes.spm.dir);
SPM.load;

EL = eeglabClass(aap.directory_conventions.toolboxes.eeglab.dir,'requiredPlugins',strsplit(aap.directory_conventions.toolboxes.eeglab.extraparameters.requiredPlugins,':'));
EL.load;
CHANNELFILE = fullfile(EL.dipfitPath,'standard_BESA','standard-10-5-cap385.elp');
EL.close;

% SITE-SPECIFIC CONFIGURATION:
aap.options.wheretoprocess = 'qsub_nonDCS'; % queuing system			% typical value localsingle or qsub_nonDCS
aap.options.aaparallel.numberofworkers = 4;
aap.options.aaparallel.memory = 4;
aap.options.aaworkerGUI = 0;

%% PIPELINE
% Directory & sub-directory for analysed data:
aap.acq_details.root = fullfile(aap.acq_details.root,'aa_demo');
aap.directory_conventions.analysisid = 'lemon'; 

% Pipeline customisation
aap = aas_addinitialstream(aap,'channellayout',{CHANNELFILE});

aap.tasksettings.aamod_meeg_converttoeeglab.removechannel = 'VEOG';
aap.tasksettings.aamod_meeg_converttoeeglab.downsample = 250;
aap.tasksettings.aamod_meeg_converttoeeglab.diagnostics.freqrange = [1 120];
aap.tasksettings.aamod_meeg_converttoeeglab.diagnostics.freq = [6 10 50];

aap.tasksettings.aamod_meeg_filter.hpfreq = 1;
aap.tasksettings.aamod_meeg_filter.bsfreq = cell2mat(arrayfun(@(x) [x-5 x+5]', [50 100], 'UniformOutput', false))';
aap.tasksettings.aamod_meeg_filter.diagnostics = aap.tasksettings.aamod_meeg_converttoeeglab.diagnostics;

aap.tasksettings.aamod_meeg_cleanartifacts.criteria.Highpass = 'off';
aap.tasksettings.aamod_meeg_cleanartifacts.criteria.LineNoiseCriterion = 'off';
aap.tasksettings.aamod_meeg_cleanartifacts.criteria.FlatlineCriterion = 5; % maximum tolerated flatline duration in seconds
aap.tasksettings.aamod_meeg_cleanartifacts.criteria.ChannelCriterion = 0.8; % minimum channel correlation
aap.tasksettings.aamod_meeg_cleanartifacts.criteria.BurstCriterion = 20; % 5 (recommended by Makoto's pres) is too agressive; 10 to *20* (according to the evaluation paper)
aap.tasksettings.aamod_meeg_cleanartifacts.criteria.Distance = 'riemannian'; % Riemann adapted processing is a newer method to estimate covariance matrices
aap.tasksettings.aamod_meeg_cleanartifacts.criteria.BurstRejection = 'off'; % correcting data using ASR instead of removing
aap.tasksettings.aamod_meeg_cleanartifacts.criteria.WindowCriterion = 0.25; % if more than this % of channels still show above-threshold amplitudes, reject this window (0.05 - 0.3)
aap.tasksettings.aamod_meeg_cleanartifacts.interpolate = 'spherical';

aap.tasksettings.aamod_meeg_rereference.reference = 'average';
aap.tasksettings.aamod_meeg_rereference.diagnostics = aap.tasksettings.aamod_meeg_converttoeeglab.diagnostics;

aap.tasksettings.aamod_meeg_ica.PCA = 'rank';
aap.tasksettings.aamod_meeg_ica.iterations = 2000;
aap.tasksettings.aamod_meeg_ica.method = 'AMICA';
aap.tasksettings.aamod_meeg_ica.options.AMICA.num_models = 1; % learn 1 model
% reject outliers (>3 SD) for the first 15 iterations 
aap.tasksettings.aamod_meeg_ica.options.AMICA.numrej = 15; 
aap.tasksettings.aamod_meeg_ica.options.AMICA.rejint = 1;
aap.tasksettings.aamod_meeg_ica.options.AMICA.rejsig = 3;

aap.tasksettings.aamod_meeg_dipfit.transformation = CHANNELFILE;
aap.tasksettings.aamod_meeg_dipfit.volumeCondutionModel = fullfile('standard_BESA','standard_BESA.mat');
aap.tasksettings.aamod_meeg_dipfit.mri = fullfile('standard_BEM','standard_mri.mat');
aap.tasksettings.aamod_meeg_dipfit.rejectionThreshold = 100; % keep all
aap.tasksettings.aamod_meeg_dipfit.constrainSymmetrical = 1;

% Automatic IC rejection using ICLabel label probability (brain > 0.7) and and residual variance (< 0.15) from dipole fitting (if performed).
aap.tasksettings.aamod_meeg_icclassification.method = 'ICLabel';
aap.tasksettings.aamod_meeg_icclassification.criteria.probBrain = 0.7;
aap.tasksettings.aamod_meeg_icclassification.criteria.rv = 0.15;

%% DATA
% Directory for raw data:
aap.directory_conventions.meegsubjectoutputformat = 'sub-%06d';
aap.directory_conventions.subject_directory_format = 1;
aarawdir = strsplit(aap.directory_conventions.rawmeegdatadir,':');
aap.directory_conventions.rawmeegdatadir = fullfile(aarawdir{cell_index(aarawdir,'aa_demo')},'lemon');
aas_makedir(aap,aap.directory_conventions.rawmeegdatadir);
for subj = SUBJS
    if exist(fullfile(aap.directory_conventions.rawmeegdatadir,sprintf(aap.directory_conventions.meegsubjectoutputformat,subj)),'dir'), continue; end
    aas_shell(['cd ' aap.directory_conventions.rawmeegdatadir ';wget -r --no-host-directories --directory-prefix=' sprintf(aap.directory_conventions.meegsubjectoutputformat,subj) ' --cut-dirs=7 ftp://ftp.gwdg.de/pub/misc/MPI-Leipzig_Mind-Brain-Body-LEMON/EEG_MPILMBB_LEMON/EEG_Raw_BIDS_ID/' sprintf(aap.directory_conventions.meegsubjectoutputformat,subj)]);
end

% Add subject (full):
aap = aas_add_meeg_session(aap,'run1');
for subj = SUBJS
    eegacq = cellstr(spm_file(spm_select('FPListRec',meeg_findvol(aap,subj,'fullpath',true),'.*vhdr'),'filename'));
    if numel(eegacq) ~= numel(aap.acq_details.meeg_sessions), aas_log(aap,false,'The numbers of EEG sessions and EEG acquisitions do not match'); end
    aap = aas_addsubject(aap,{subj []},'functional',eegacq);
end

%% RUN
aa_doprocessing(aap);
% aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));