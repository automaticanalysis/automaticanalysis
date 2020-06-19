% Automatic analysis
% User master example (aa version 5.*.*)
%
% This user script demonstrates a basic diffusion ROI pipeline with diffusion kurtosis
% imaging (DKI).
%
% For internal use at MRC CBU, Cambridge, UK - requires access to the CBU imaging
% system.
%
% v2: Johan Carlin, MRC CBU, 07-08-2018
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
aap=aarecipe('aap_tasklist_diffusion.xml');

% this example uses SPM tools in the user script, so we have to ensure SPM is
% on the path
spmhit = which('spm_spm');
if any(spmhit)
    assert(strcmp(fileparts(spmhit), aap.directory_conventions.toolboxes.spm.dir), ...
        'spm on path differs from aap.directory_conventions.toolboxes.spm.dir');
else
    fprintf('adding spmdir to path: %s\n', aap.directory_conventions.toolboxes.spm.dir);
    SPM = spmClass(aap.directory_conventions.toolboxes.spm.dir);
    SPM.load;
end

aap = aas_renamestream(aap,'aamod_diffusion_dartel_denormDKI_00001','grey','normalised_white', 'input');
aap = aas_renamestream(aap,'aamod_diffusion_dartel_denormDKI_00001','grey','native_white','output');

% Modify standard recipe module selection here if you'd like
aap.options.wheretoprocess = 'qsub'; % queuing system	% typical value localsingle or qsub

aap.options.autoidentifystructural_chooselast = 1;
aap.tasksettings.aamod_segment8.writenormimg=0;
aap.tasksettings.aamod_dartel_norm_write.fwhm=1;

aap.tasksettings.aamod_diffusion_bet.bet_f_parameter = 0.2;
aap.tasksettings.aamod_diffusion_smooth.FWHM = 2.5;        
aap.tasksettings.aamod_diffusion_dartel_denormDKI.interp=4;

%% STUDY
aap=aas_addinitialstream(aap,'rois',{...
    '/imaging/local/software/AA/test_resources/diffusion/craddock_ROI_841_Linda_FCpaper.nii'});

% Directory for analysed data - appending to cbu default root directory
aap.acq_details.root = fullfile(aap.acq_details.root,'aa_demo');
aap.directory_conventions.analysisid = 'diffusion'; 
aap.directory_conventions.subject_directory_format = 3;

% Add data
aap = aas_add_diffusion_session(aap,'diffusion');
for s = 1:size(SUBJ,1)
   diffser = sscanf(basename(spm_select('FPListRec',mri_findvol(aap,SUBJ{s,2},1),'dir','.*_CBU_DKI_30dir_2bvals$')),aap.directory_conventions.seriesoutputformat);
   aap = aas_addsubject(aap,SUBJ{s,1},SUBJ{s,2},'diffusion',diffser);
end

%% DO ANALYSIS
aa_doprocessing(aap);
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
