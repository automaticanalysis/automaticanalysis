function [aap,resp]=aamod_freesurfer_initialise(aap,task,subj)

resp='';

switch task
    case 'report'
        
    case 'doit'
        
        % Load the Structural
        Simg = aas_getfiles_bystream(aap,subj,'structural');
        if size(Simg,1) > 1
            aas_log(aap, false, sprintf('Found more than 1 structural images, using structural %d', ...
                aap.tasklist.currenttask.settings.structural));
            Simg = Simg(aap.tasklist.currenttask.settings.structural, :);
        end
        Simg = strvcat2cell(Simg);
        
        % Set subject paths
        subjname = aap.acq_details.subjects(subj).subjname;
        subjpath = aas_getsubjpath(aap,subj);
        
        setenv('SUBJECTS_DIR', fileparts(subjpath))
        setenv('FREESURFER_DIR', aap.directory_conventions.freesurferdir)
        
        %% initialise fileserver folder structure and nii and mgh files
        aas_freesurfer_init(aap, subjpath, subjname, Simg, 1);
        
        %% Try to delete old freesurfer running flags
        if exist(fullfile(subjpath, 'ANAT', subjname, 'scripts', 'IsRunning.lh+rh'), 'file')
            unix(['rm ' fullfile(subjpath, 'ANAT', subjname, 'scripts', 'IsRunning.lh+rh')]);
        end
        
        %%  make output stream
        % (JC: now specific to FS directories rather than everything under
        % the sun)
        subdir = aas_getsubjpath(aap,subj);
        outs = [];
        for d = {'RAW','ANAT','mri'}
            outs = [outs; dirrec(fullfile(subdir,d{1}))];
        end
        aap = aas_desc_outputs(aap,subj,'freesurfer',outs);
end
