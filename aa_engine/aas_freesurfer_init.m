function aas_freesurfer_init(aap, subj_dir, subj_id, Simg, OverwriteFiles)
%INIT Summary of this function goes here
%   Detailed explanation goes here
% subj_dir is directory of freesurfer processing
% subj_id is the subject ID, used by freesurfer
% Simg is the structural image used by freesurfer
% OverwriteFiles = 1-->when finding nii, mgh, enz. files they will be recreated, if 0
% will leave them alone

global IDs;
IDs = 1;

% Nifti dir
if ~exist(fullfile(subj_dir, 'RAW', subj_id, 'NII'), 'dir')
    mkdir(fullfile(subj_dir, 'RAW', subj_id, 'NII'));
end

% Logs dir
if ~exist(fullfile(subj_dir, 'LOGS', subj_id), 'dir')
    mkdir(fullfile(subj_dir, 'LOGS', subj_id));
end

%% Make FreeSurfer directory structure and plant converted anatomicals

if ~exist(fullfile(subj_dir, 'ANAT', subj_id),'dir')
    mkdir(fullfile(subj_dir, 'ANAT', subj_id));
end

SkipAnatConvert = false;

if ~isempty(dir([subj_dir, 'ANAT', subj_id, 'mri', 'orig', '*.mgz']))
    if ~OverwriteFiles
        SkipAnatConvert = true;
    end
end

if ~SkipAnatConvert
    if ~exist(fullfile(subj_dir, 'ANAT', subj_id, 'mri', 'orig'), 'dir')
        mkdir(fullfile(subj_dir, 'ANAT', subj_id, 'mri', 'orig'));
    end
    if ~exist(fullfile(subj_dir, 'mri'), 'dir')
        mkdir(fullfile(subj_dir, 'mri'));
    end
    
    for anat_idx = 1:length(Simg)
        thisanat = Simg{anat_idx};
        newanat = fullfile(subj_dir, 'ANAT', subj_id, 'mri', 'orig', sprintf('%03d.mgz' , anat_idx));
        [s,w] = aas_runFScommand(aap,['mri_convert ', thisanat, ' ', newanat]);
        newanat = fullfile(subj_dir, 'mri', sprintf('%03d.mgz' , anat_idx));
        [s,w] = aas_runFScommand(aap,['mri_convert ', thisanat, ' ', newanat]);
        fid = fopen(fullfile(subj_dir, 'LOGS', subj_id, sprintf('nii2mgh_%03d.log' , anat_idx)), 'w');
        fprintf(fid,w);
        fclose(fid);
        if s==1
            aas_log(aap,true,'ERROR: Conversion of anatomicals from NIFTI to MGH failed. Please consult the nii2mgh log for more information. Exiting.');
        end
    end
else
    aas_log(aap,false,'Conversion of anatomicals from NIFTI to MGH was skipped, existing MGHs will be used.');
end

aas_log(aap,false,['Initialization freesurfer for subject ', subj_id, ' ran succesfully!']);

end
