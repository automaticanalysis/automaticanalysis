function [aap,resp]=aamod_freesurfer_register(aap,task,subj)

resp='';

switch task
    case 'report'
        
    case 'doit'
        
        % Set subject paths
        subjpath = aas_getsubjpath(aap,subj);
        subjname = basename(subjpath);
        
        setenv('SUBJECTS_DIR', fileparts(subjpath))
        
        %% Try to delete old freesurfer running flags
        % Get data
        inputstreams = aas_getstreams(aap,'in');
        struct = aas_getfiles_bystream(aap, subj,inputstreams{2});
        reg = fullfile(subjpath,aap.directory_conventions.structdirname,'register.dat');        
        if numel(inputstreams) == 2
            srcstream{1} = inputstreams{2};
        else
            for i = 3:numel(inputstreams)
                srcstream{i-2} = inputstreams{i};
            end
        end
        FWHM = 0;
        if isfield(aap.tasklist.currenttask.settings, 'FWHM')
            FWHM = aap.tasklist.currenttask.settings.FWHM;
        end
        
        
        % Register T1 to FS
        switch inputstreams{2}
            case 'structural'
                mod = 't1';
        end
        FScommand = sprintf('bbregister --s %s --mov %s --reg %s --init-spm --%s',subjname,struct,reg,mod);
        [s, w] = aas_runFScommand(aap,FScommand);
        
        % Register SRC to FS
        flags = aap.spm.defaults.coreg.write;
        flags.interp = aap.tasklist.currenttask.settings.interp;
        flags.which = [1 0];
        srcwfs = '';
        for i = 1:numel(srcstream)
            fsrc = aas_getfiles_bystream(aap, subj,srcstream{i}); 
            for f = 1:size(fsrc,1)
                src = deblank(fsrc(f,:));
                spm_reslice({struct src},flags);
                src = fullfile(fileparts(src),[flags.prefix basename(src)]);
                FScommand = sprintf('mri_vol2surf --mov %s.nii --reg %s --hemi rh --o %s2FS_rh.mgh',src,reg,src);        
                if FWHM, FScommand = [FScommand sprintf(' --surf-fwhm %d',FWHM)]; end            
                [s, w] = aas_runFScommand(aap,FScommand);
                FScommand = sprintf('mri_vol2surf --mov %s.nii --reg %s --hemi lh --o %s2FS_lh.mgh',src,reg,src);        
                if FWHM, FScommand = [FScommand sprintf(' --surf-fwhm %d',FWHM)]; end            
                [s, w] = aas_runFScommand(aap,FScommand);
                %% Output stream
                srcwfs(end+1,:) = [strrep(src,[subjpath '/'],'') '2FS_rh.mgh'];
                srcwfs(end+1,:) = [strrep(src,[subjpath '/'],'') '2FS_lh.mgh'];
            end
            aap = aas_desc_outputs(aap,subj,[srcstream{i} '_FS'],srcwfs);
        end
end
