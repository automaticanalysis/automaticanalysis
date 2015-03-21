function [aap,resp]=aamod_freesurfer_autorecon(aap,task,subj)

resp='';

switch task
    case 'report'
        
    case 'doit'
        
        % Set subject paths
        subjpath = aas_getsubjpath(aap,subj);
        subjname = basename(subjpath);
        
        setenv('SUBJECTS_DIR', fileparts(subjpath))
        
        %% Try to delete old freesurfer running flags
        if exist(fullfile(subjpath, 'ANAT', subjname, 'scripts', 'IsRunning.lh+rh'), 'file')
            unix(['rm ' fullfile(subjpath, 'ANAT', subjname, 'scripts', 'IsRunning.lh+rh')]);
        end
        
        FScommand = ['recon-all -subjid ' subjname ' ' aap.tasklist.currenttask.settings.extraoptions];
        
        disp(FScommand)
        
        [s, w] = aas_runFScommand(aap,FScommand);
                
        if s==1 %|| ~isempty(strfind(w, 'ERROR'))
            disp(w);
            error('Some system ERROR');
        end
        
        if aap.tasklist.currenttask.settings.verbose
            disp(w);
        end
        
        %%  make output stream
        % now specific to freesurfer dirs. 
        subdir = fullfile(aas_getsubjpath(aap, subj)); % freesurfer autorecon1 dir
        outs = [];
        dirs = {'ANAT','bem','label','mri','scripts','src','stats',...
            'surf','touch',};
        for d = dirs
            outs = [outs dirrec(fullfile(subdir,d{1}))];
        end
        aap = aas_desc_outputs(aap,subj,'freesurfer',outs);
end
