% Automatic analysis
% User master script example (aa version 5.*.*)
%
% This is an example how to process BIDS multimodal NIfTI dataset "7t_trt" 
% (https://github.com/INCF/BIDS-examples/tree/master/7t_trt
% It uses a tasklist BIDS7T_tasklist.xml (included)
%
% Tibor Auer, MRC-CBSU
% 08-02-2016

%% INITIALISE
clear

aa_ver5

%% DEFINE SPECIFIC PARAMETERS
%  Default recipe without model
aap=aarecipe('aap_parameters_defaults_CBSU.xml','BIDS7T_tasklist.xml');
aap = aas_configforSPM12(aap);

% Modify standard recipe module selection here if you'd like
aap.options.wheretoprocess = 'qsub'; % queuing system			% typical value qsub | localsingle
aap.options.NIFTI4D = 1;										% typical value 0
aap.options.email='xy00@mrc-cbu.cam.ac.uk';

aap.tasksettings.aamod_dartel_norm_write.vox = 1;
for b = 1:4
    aap.tasksettings.aamod_slicetiming(b).autodetectSO = 1;
    aap.tasksettings.aamod_slicetiming(b).refslice = 16;
    aap.tasksettings.aamod_norm_write_dartel(b).vox = [3 3 3];
    aap.tasksettings.aamod_smooth(b).FWHM = 5;
end

%% STUDY
% Directory for analysed data
aap.acq_details.root = '/imaging/xy00/aa'; 
aap.directory_conventions.analysisid = 'BIDS_7T'; 

% Add data
aap.directory_conventions.rawdatadir = '/Data/BIDS/7t_trt';
aap.acq_details.numdummies = 0;
aap.acq_details.input.combinemultiple = 1;
aap.options.autoidentifystructural_choosefirst = 1;
aap = aas_processBIDS(aap);

%% DO ANALYSIS
aa_doprocessing(aap);
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));