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
        
        aas_log(aap,false,FScommand)
        
        [s, w] = aas_runFScommand(aap,FScommand,~aap.tasklist.currenttask.settings.verbose);
                
        if s==1 %|| ~isempty(strfind(w, 'ERROR'))
            aas_log(aap,true,w);
        end
        
        %%  make output stream
        % now specific to freesurfer dirs. 
        subjpath = fullfile(aas_getsubjpath(aap, subj)); % freesurfer autorecon1 dir
        outs = cellfun(@(x) spm_select('FPListRec',fullfile(subjpath,x),'.*'),...
            {'ANAT','bem','label','mri','scripts','src','stats', 'surf','touch',},'UniformOutput',false);
        outs = char(outs(cellfun(@(x) ~isempty(x), outs)));
        aap = aas_desc_outputs(aap,subj,'freesurfer',outs);
end
