function [aap,resp]=aamod_freesurfer_register(aap,task,subj)

resp='';

switch task
    case 'report'
        
    case 'doit'
        % Settings
        srcstream = aas_getstreams(aap,'input');
        FWHM = 0;
        if isfield(aap.tasklist.currenttask.settings, 'FWHM')
            FWHM = aap.tasklist.currenttask.settings.FWHM;
        end        
        resliceflags = aap.spm.defaults.coreg.write;
        resliceflags.interp = aap.tasklist.currenttask.settings.interp;
        resliceflags.which = [1 0];
        
        % Set paths and data
        switch aap.tasklist.currenttask.domain
            case 'study'
                localroot = aas_getstudypath(aap);
                setenv('SUBJECTS_DIR', localroot)
                system(sprintf('ln -s %s/subjects/fsaverage %s/fsaverage',aap.directory_conventions.freesurferdir,localroot));
                system(sprintf('ln -s %s/subjects/lh.EC_average %s/lh.EC_average',aap.directory_conventions.freesurferdir,localroot));
                system(sprintf('ln -s %s/subjects/rh.EC_average %s/rh.EC_average',aap.directory_conventions.freesurferdir,localroot));
                subjname = 'fsaverage';
                aas_makedir(aap,fullfile(localroot,aap.directory_conventions.structdirname));
                % get template try FSL's
                struct = fullfile(aap.directory_conventions.fsldir,'data','standard','MNI152_T1_1mm.nii.gz'); % use FSL highres
                if exist(struct,'file')
                    gunzip(struct,fullfile(localroot,aap.directory_conventions.structdirname)); 
                    struct = fullfile(fullfile(localroot,aap.directory_conventions.structdirname),'MNI152_T1_1mm.nii');
                else % use SPM's
                    struct = fullfile(aap.directory_conventions.spmdir,aap.directory_conventions.T1template);
                    copyfile(struct,fullfile(localroot,aap.directory_conventions.structdirname));
                    struct = spm_file(struct,'path',fullfile(localroot,aap.directory_conventions.structdirname));
                end
                esrc = aas_getfiles_bystream(aap,'study', [], srcstream{1});
                spm_reslice({esrc(1,:) struct},resliceflags);
                struct = spm_file(struct,'prefix',resliceflags.prefix);
                mod = 't1';
                indices = [];
            case 'subject'
                localroot = aas_getsubjpath(aap,subj);
                setenv('SUBJECTS_DIR', fileparts(localroot))
                subjname = basename(localroot);
                struct = aas_getfiles_bystream(aap,subj,srcstream{2});
                srcstr = strsplit(srcstream{2},'.');
                switch srcstr{end}
                    case 'structural'
                        mod = 't1';
                end
                srcstream(1:2) = []; % keep sources only
                indices = subj;
        end
        
        reg = fullfile(localroot,aap.directory_conventions.structdirname,'register.dat');        
        
        % Register T1 to FS        
        FScommand = sprintf('bbregister --s %s --mov %s --reg %s --init-spm --spm-nii --%s',subjname,struct,reg,mod);
        [s, w] = aas_runFScommand(aap,FScommand);
        
        % Register SRC to FS
        resliceflags = aap.spm.defaults.coreg.write;
        resliceflags.interp = aap.tasklist.currenttask.settings.interp;
        resliceflags.which = [1 0];
        srcwfs = {};
        for i = 1:numel(srcstream)
            fsrc = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain, indices, srcstream{i}); 
            for f = 1:size(fsrc,1)
                src = deblank(fsrc(f,:));
                spm_reslice({struct src},resliceflags);
                src = fullfile(fileparts(src),[resliceflags.prefix basename(src)]);
                FScommand = sprintf('mri_vol2surf --mov %s.nii --reg %s --hemi rh --o %s2FS_rh.mgh',src,reg,src);        
                if FWHM, FScommand = [FScommand sprintf(' --surf-fwhm %d',FWHM)]; end            
                [s, w] = aas_runFScommand(aap,FScommand);
                FScommand = sprintf('mri_vol2surf --mov %s.nii --reg %s --hemi lh --o %s2FS_lh.mgh',src,reg,src);        
                if FWHM, FScommand = [FScommand sprintf(' --surf-fwhm %d',FWHM)]; end            
                [s, w] = aas_runFScommand(aap,FScommand);
                %% Output stream
                srcwfs{end+1} = [strrep(src,[localroot '/'],'') '2FS_rh.mgh'];
                srcwfs{end+1} = [strrep(src,[localroot '/'],'') '2FS_lh.mgh'];
            end
            aap = aas_desc_outputs(aap,aap.tasklist.currenttask.domain,indices,[srcstream{i} '_FS'],char(srcwfs));
        end
end
