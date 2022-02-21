function aatest_lemon_meegsource(parameterfile, deleteprevious, wheretoprocess)
% This script tests an advanced EEG pipeline on a single participant in the LEMON dataset: http://fcon_1000.projects.nitrc.org/indi/retro/MPI_LEMON.html.
% The retrieval of the dataset requires aws-cli
%
% It requires the following software to be configured in your parameterset
%   - SPM
%   - FSL
%   - EEGLAB with extensions Fileio, bva-io, clean_rawdata, AMICA, dipfit, Fieldtrip-lite, firfilt, fitTwoDipoles, ICLabel, Viewprops
%   - FieldTrip
% See aa_parametersets/aap_parameters_defaults_UoS.xml, lines 21-56 for example configuration

aap = aa_test_inittest(mfilename('fullpath'), parameterfile, deleteprevious, wheretoprocess);

SUBJS = [ 32301 ];

%% RECIPE
EL = aas_inittoolbox(aap,'eeglab');
EL.load;
CHANNELFILE = fullfile(EL.dipfitPath,'standard_BESA','standard-10-5-cap385.elp');
EL.close;

MRIFILE = fullfile(aap.directory_conventions.fsldir,'data/standard/MNI152_T1_1mm.nii.gz');

%% PIPELINE
% Directory & sub-directory for analysed data:
aap.acq_details.root = fullfile(aap.acq_details.root,'aa_demo');
aap.directory_conventions.analysisid = 'lemon'; 

% Pipeline customisation
aap = aas_addinitialstream(aap,'channellayout',{CHANNELFILE});
aap = aas_addinitialstream(aap,'MNI_1mm',{MRIFILE});
aap.tasksettings.aamod_importfilesasstream(2).unzip = 'gunzip';

aap.tasksettings.aamod_structuralfromnifti.sfxformodality = 'T1w'; % suffix for structural
aap.tasksettings.aamod_segment8.combine = [0.05 0.05 0.05 0.05 0.5 0];
aap.tasksettings.aamod_segment8.writenormimg = 0; % write normialised structural
aap = aas_renamestream(aap,'aamod_coreg_general_00001','reference','MNI_1mm','input');
aap = aas_renamestream(aap,'aamod_coreg_general_00001','input','structural','input');
aap = aas_renamestream(aap,'aamod_coreg_general_00001','output','structural','output');
aap.tasksettings.aamod_meeg_prepareheadmodel.method = 'simbio';
aap.tasksettings.aamod_meeg_prepareheadmodel.options.simbio.downsample = 2;
aap.tasksettings.aamod_meeg_prepareheadmodel.options.simbio.meshshift = 0.1;
aap.tasksettings.aamod_meeg_preparesourcemodel.method = 'grid';
aap.tasksettings.aamod_meeg_preparesourcemodel.options.grid.resolution = '10';
aap = aas_renamestream(aap,'aamod_norm_write_00001','structural','aamod_coreg_general_00001.structural','input');
aap = aas_renamestream(aap,'aamod_norm_write_00001','epi','aamod_coreg_general_00001.structural','input');
aap = aas_renamestream(aap,'aamod_norm_write_00001','epi','structural','output');
aap.tasksettings.aamod_norm_write.bb = [-90 90 -126 91 -72 109];
aap.tasksettings.aamod_norm_write.vox = [1 1 1];

aap.tasksettings.aamod_meeg_converttoeeglab.removechannel = 'VEOG';
aap.tasksettings.aamod_meeg_converttoeeglab.downsample = 250;
aap.tasksettings.aamod_meeg_converttoeeglab.diagnostics.freqrange = [1 120];
aap.tasksettings.aamod_meeg_converttoeeglab.diagnostics.freq = [6 10 50];
% - correct/harmonise events
aap.tasksettings.aamod_meeg_converttoeeglab.toEdit(1).subject = '*';
aap.tasksettings.aamod_meeg_converttoeeglab.toEdit(1).session = '*';
aap.tasksettings.aamod_meeg_converttoeeglab.toEdit(1).event(1) = struct('type','^S.*','operation','keep'); % keep only events starting with 'S'
aap.tasksettings.aamod_meeg_converttoeeglab.toEdit(1).event(2) = struct('type','S  1','operation','unique:last'); % remove duplicates of 'S  1'
aap.tasksettings.aamod_meeg_converttoeeglab.toEdit(1).event(3) = struct('type','S  1','operation','iterate'); % 'S  1' -> 'S  101', 'S  102', 'S  103',...
aap.tasksettings.aamod_meeg_converttoeeglab.toEdit(2).subject = '*';
aap.tasksettings.aamod_meeg_converttoeeglab.toEdit(2).session = '*';
aap.tasksettings.aamod_meeg_converttoeeglab.toEdit(2).event(1) = struct('type','S  101','operation','ignorebefore'); % remove heading samples

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
aap.tasksettings.aamod_meeg_dipfit.rejectionThreshold = 100; % keep all
aap.tasksettings.aamod_meeg_dipfit.constrainSymmetrical = 1;

% Automatic IC rejection using ICLabel label probability (brain > 0.7) and and residual variance (< 0.15) from dipole fitting (if performed).
aap.tasksettings.aamod_meeg_icclassification.method = 'ICLabel';
aap.tasksettings.aamod_meeg_icclassification.criteria.prob = 'Brain>0.7'; % Eye<0.8:*Muscle<0.8
aap.tasksettings.aamod_meeg_icclassification.criteria.rv = 0.15;

for b = 1:2
    aap.tasksettings.aamod_meeg_timefrequencyanalysis(b).ignorebefore = -6; % ignore the first 5 trials for each trialtype 
    aap.tasksettings.aamod_meeg_timefrequencyanalysis(b).timefrequencyanalysis.method = 'mtmfft';
    aap.tasksettings.aamod_meeg_timefrequencyanalysis(b).timefrequencyanalysis.taper = 'hanning';
    aap.tasksettings.aamod_meeg_timefrequencyanalysis(b).timefrequencyanalysis.foi = [1:0.5:13 15 20 25 32 40 60 70 80 95 110 120];
    aap.tasksettings.aamod_meeg_timefrequencyanalysis(b).diagnostics.snapshotfwoi = [...
        1 3.5;... % delta
        4 7.5;... % theta
        8 13;... % alpha
        14 32;... % beta
        33 80;... % low-gamma
        81 120;... % high-gamma
        ];
    aap.tasksettings.aamod_meeg_timefrequencyanalysis(b).contrastoperation = 'ratio';
end
aap.tasksettings.aamod_meeg_timefrequencyanalysis(1).weightedaveraging = 1;
aap.tasksettings.aamod_meeg_timefrequencyanalysis(2).diagnostics.snapshottwoi = [[0:120000:7*120000]' [0:120000:7*120000]'+120000];

for b = 1:2
    aap = aas_renamestream(aap,sprintf('aamod_meeg_sourcereconstruction_%05d',b),'input','timefreq');
    aap.tasksettings.aamod_meeg_sourcereconstruction(b).realignelectrodes.target = 'scalp';
    aap.tasksettings.aamod_meeg_sourcereconstruction(b).realignelectrodes.method = 'spherefit';
    aap.tasksettings.aamod_meeg_sourcereconstruction(b).diagnostics = struct_update(aap.tasksettings.aamod_meeg_sourcereconstruction(b).diagnostics,...
        aap.tasksettings.aamod_meeg_timefrequencyanalysis(b).diagnostics,'Mode','update');
end

%% DATA
% Directory for raw data:
aap.directory_conventions.subject_directory_format = 1;
% - eeg
aap.directory_conventions.meegsubjectoutputformat = 'sub-%06d/RSEEG';
% - mri
aap.directory_conventions.subjectoutputformat = 'sub-%06d';

aap = aas_add_meeg_session(aap,'run1');
for subj = SUBJS
    eegacq = cellstr(spm_file(spm_select('FPListRec',meeg_findvol(aap,subj,'fullpath',true),'.*vhdr'),'filename')); % filename only for EEG
    mriacq = cellstr(spm_select('FPListRec',mri_findvol(aap,subj,'fullpath',true),'.*_ses-01_acq-mp2rage_T1w.nii.gz')); % fullfile for MRI
    aap = aas_addsubject(aap,{subj subj},'structural',mriacq,'functional',eegacq);
end

%% Epoching
aap = aas_add_meeg_event(aap,'aamod_meeg_epochs','*','run1','segment-1','S  101:S  103',0);
aap = aas_add_meeg_event(aap,'aamod_meeg_epochs','*','run1','segment-2','S  103:S  105',0);
aap = aas_add_meeg_event(aap,'aamod_meeg_epochs','*','run1','segment-3','S  105:S  107',0);
aap = aas_add_meeg_event(aap,'aamod_meeg_epochs','*','run1','segment-4','S  107:S  109',0);
aap = aas_add_meeg_event(aap,'aamod_meeg_epochs','*','run1','segment-5','S  109:S  111',0);
aap = aas_add_meeg_event(aap,'aamod_meeg_epochs','*','run1','segment-6','S  111:S  113',0);
aap = aas_add_meeg_event(aap,'aamod_meeg_epochs','*','run1','segment-7','S  113:S  115',0);
aap = aas_add_meeg_event(aap,'aamod_meeg_epochs','*','run1','segment-8','S  115:end',0);

aap = aas_add_meeg_event(aap,'aamod_meeg_epochs','*','run1','EC','S210',0);
aap = aas_add_meeg_event(aap,'aamod_meeg_epochs','*','run1','EO','S200',0);
aap.tasksettings.aamod_meeg_epochs.timewindow = [0 2000];

%% Analysis
aap = aas_add_meeg_trialmodel(aap,'aamod_meeg_timefrequencyanalysis_00001','*','singlesession:run1','+1xEC','avg','ECAVG');
aap = aas_add_meeg_trialmodel(aap,'aamod_meeg_timefrequencyanalysis_00001','*','singlesession:run1','+1xEO','avg','EOAVG');
aap = aas_add_meeg_trialmodel(aap,'aamod_meeg_timefrequencyanalysis_00001','*','singlesession:run1','+1xEC|-1xEO','avg','ECAVGminusEOAVG');

aap = aas_add_meeg_trialmodel(aap,'aamod_meeg_timefrequencyanalysis_00002','*','singlesession:run1','+1xEC','segmentavg','ECCONT');
aap = aas_add_meeg_trialmodel(aap,'aamod_meeg_timefrequencyanalysis_00002','*','singlesession:run1','+1xEO','segmentavg','EOCONT');

%% RUN
aa_doprocessing(aap);
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));