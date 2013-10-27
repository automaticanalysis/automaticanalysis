% Automatic analysis user script

addpath('/system/cn_developer/camneuro/release-beta-0.0/automaticanalysis/aa_engine')
addpath('/imaging/rhodri/hAAckathon/dti_scripts');
addpath('/imaging/leire/DTI_CamCan/Scripts');

% Does two iterations, as at present, bedpostx won't run when launched from
% within a condor module

% for iteration=iteration_start:iteration_end
%      if iteration==1
%         aap=aarecipe('aap_tasklist_dti_part2_Leire_withprobtracx_3_TEST.xml');
%         aap.options.wheretoprocess='localsingle';
%     elseif iteration==2
        aap=aarecipe('aap_tasklist_diffusion.xml');
         if exist('runlocally','var') && runlocally
             aap.options.wheretoprocess='localsingle';
         else
             aap.options.wheretoprocess='condor';
         end;
%     end;
    

% DEFINE STUDY SPECIFIC PARAMETERS
aap.options.aa_minver=4.0;
aap.directory_conventions.subject_directory_format=3;



% Where to put the analyzed data
aap.acq_details.root = '/imaging/leire/DTI_CamCan';
aap.directory_conventions.analysisid='dti_analysis';
aap.options.autoidentifyfieldmaps=false;
aap.directory_conventions.rawdatadir='/imaging/rhodri/camcan/interim/raw/raw/cc700-interim-raw/raw';

% aap=aas_addsubject(aap,'103515'); % files are in .nii not .nii.gz
% format
fn=dir('/imaging/rhodri/camcan/interim/raw/raw/cc700-interim-raw/raw');
fnind=1;
for fnind=1:length(fn)
    if (length(aap.acq_details.subjects)==1)
        break
    end;
    if (length(fn(fnind).name)==6) && exist(fullfile('/imaging/leire/DTI_CamCan',fn(fnind).name),'dir')
        aap=aas_addsubject(aap,fn(fnind).name);
    end;
end;
aap=aas_addsubject(aap,'CC110037_CBU110544',[],[],[],16);

% One DTI session
aap=aas_add_diffusion_session(aap,'diffusion');

% PICK ONE SUBJECT FOR TESTING
%aap.acq_details.subjects=aap.acq_details.subjects(1);

% DO PROCESSING
aa_doprocessing(aap);

% end;
