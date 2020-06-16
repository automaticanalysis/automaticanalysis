% Automatic analysis
% User master example (aa version 5.*.*)
%
% This user script demonstrates a diffusion ROI pipeline with diffusion kurtosis imaging
% (DKI) and distortion correction with FSL topup.
%
% For internal use at MRC CBU, Cambridge, UK - requires access to the CBU imaging
% system.
%
% v2: Johan Carlin, MRC CBU, 08-08-2018
% v1: Tibor Auer, MRC-CBSU, 16-08-2017

%% INITIALISE
clear
aa_ver5



%% DEFINE SPECIFIC PARAMETERS
SUBJ = {...
     'S01' 160018; ...
     };

%  Default recipe without model
aap=aarecipe('aap_tasklist_diffusion_with_topup.xml');

% this example uses SPM tools in the user script, so we have to ensure SPM is
% on the path
spmhit = which('spm_spm');
if any(spmhit)
    assert(strcmp(fileparts(spmhit), aap.directory_conventions.toolboxes.spm.dir), ...
        'spm on path differs from aap.directory_conventions.toolboxes.spm.dir');
else
    fprintf('adding spmdir to path: %s\n', aap.directory_conventions.toolboxes.spm.dir);
    addpath(aap.directory_conventions.toolboxes.spm.dir);
end

aap = aas_renamestream(aap,'aamod_diffusion_dartel_denormDKI_00001','grey','normalised_white', 'input');
aap = aas_renamestream(aap,'aamod_diffusion_dartel_denormDKI_00001','grey','native_white','output');

% Modify standard recipe module selection here if you'd like
aap.options.wheretoprocess = 'qsub'; % queuing system	% typical value localsingle or qsub

aap.options.autoidentifystructural_choosefirst = 1;
aap.tasksettings.aamod_segment8.writenormimg=0;
aap.tasksettings.aamod_dartel_norm_write.fwhm=1;

aap.tasksettings.aamod_diffusion_extractnodif.tol_b0 = 10;
aap.tasksettings.aamod_diffusion_bet.bet_f_parameter = 0.2;
aap.tasksettings.aamod_diffusion_smooth.FWHM = 2.5;        
aap.tasksettings.aamod_diffusion_dartel_denormDKI.interp=4;

%% STUDY
aap=aas_addinitialstream(aap,'rois',{...
    '/imaging/local/software/AA/test_resources/diffusion/craddock_ROI_841_Linda_FCpaper.nii'});

% Directory for analysed data - appending to cbu default root directory
aap.acq_details.root = fullfile(aap.acq_details.root,'aa_demo');
aap.directory_conventions.analysisid = 'diffusion_topup'; 
aap.directory_conventions.subject_directory_format = 3;

% Add data
aap = aas_add_diffusion_session(aap,'diffusion');
for s = 1:size(SUBJ,1)
    diffser = cellfun(@(x) sscanf(x,aap.directory_conventions.seriesoutputformat),...
        cellstr(basename(spm_select('FPListRec',mri_findvol(aap,SUBJ{s,2},1),'dir','.*_diff_mbep2d_.*'))))';
    aap = aas_addsubject(aap,SUBJ{s,1},SUBJ{s,2},'diffusion',{diffser});
end

%% DO ANALYSIS
aa_doprocessing(aap);
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
