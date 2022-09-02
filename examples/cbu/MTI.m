% Automatic analysis
% User master script example (aa version 5.*.*)
%
% This script demonstrates a magnetisation transfer imaging (MTI) pipeline using the
% CamCAN dataset.
%
% For internal use at MRC CBU, Cambridge, UK - requires access to the CBU imaging
% system.
%
% v2: Johan Carlin, MRC CBU, 08-08-2018
% v1: Tibor Auer, MRC-CBSU, 08-02-2016

%% INITIALISE
clear
aa_ver5

%% DEFINE SPECIFIC PARAMETERS
SUBJ = {...
     'S01' 140905; ...
     'S02' 140910; ...
     'S03' 140913; ...
     'S04' 140928; ...
     'S05' 140931; ...
     };

aap = aarecipe('MTI.xml');

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
aap.options.wheretoprocess = 'qsub'; %'localsingle'; 

aap.options.autoidentifystructural_chooselast = 1;
aap.tasksettings.aamod_segment8.writenormimg=0;
aap.tasksettings.aamod_dartel_norm_write.fwhm=1;
aap.tasksettings.aamod_convert_specialseries.NIFTI4D = 1;

%% STUDY
% Directory for analysed data 
% (expanding off existing default root, so /imaging/$USER/aa/MTI if you haven't changed
% your defaults.)
aap.acq_details.root = fullfile(aap.acq_details.root,'aa_demo');
aap.directory_conventions.analysisid = 'mti'; 

% Add data
aap=aas_add_special_session(aap,'MTI');
for s = 1:size(SUBJ,1)
   MTser = sscanf(basename(spm_select('FPListRec',mri_findvol(aap,SUBJ{s,2},1),'dir','.*_CBU_MTR_TR50_MT$')),aap.directory_conventions.seriesoutputformat);
   aap = aas_addsubject(aap,SUBJ{s,1},SUBJ{s,2},'specialseries',{[MTser MTser+1]});
   aap = aas_addinitialstream(aap,'rois', SUBJ{s,1}, ...
     {'/imaging/camcan/templates/Juelich-maxprob-thr25-2mm.nii'});
end

%% DO ANALYSIS
aa_doprocessing(aap);
% prepare report (this can be slow...)
aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
