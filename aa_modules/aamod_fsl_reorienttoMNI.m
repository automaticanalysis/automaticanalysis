function [aap, resp] = aamod_fsl_reorienttoMNI(aap,task,subj)

resp = '';

switch task
    case 'doit'
        %% Init
        localroot = fullfile(aas_getsubjpath(aap,subj),aap.directory_conventions.structdirname);
        
        % Get the T1 template
        sTimg = aap.directory_conventions.T1template;
        if ~exist(sTimg,'file') % try in SPM
            if sTimg(1) ~= '/', sTimg = fullfile(aap.directory_conventions.spmdir,sTimg); end
        else
            sTimg = which(sTimg);
        end
        if ~exist(sTimg,'file'),
            aas_log(aap, true, sprintf('Couldn''t find template T1 image %s.', sTimg));
        end  
        
        % Get structural
        Simg = aas_getfiles_bystream(aap,subj,'structural');        
        
        %% Run
        % Obtain image geometry
        junk = spm_vol(sTimg);
        dim = abs(junk.dim);
        [junk, vox] = spm_get_bbox(sTimg); vox = abs(vox);
        bbox = dim.*vox;
        
        junk = spm_vol(Simg);
        dim = abs(junk.dim);
        [junk, vox] = spm_get_bbox(Simg); vox = abs(vox);
        dim(1) = ceil(bbox(1)/vox(1));
        dim(2) = ceil(bbox(2)/vox(2));
        dim(3) = ceil(bbox(3)/vox(3));
        
        % Reslice template to image geometry
        aas_runfslcommand(aap,sprintf('fslcreatehd %d %d %d 1 %8.6f %8.6f %8.6f 1 0 0 0 2 %s',dim(:),vox(:),fullfile(localroot,'tmp.nii')));
        aas_runfslcommand(aap,sprintf('flirt -in %s -applyxfm -init %s -out %s -paddingsize 0.0 -interp trilinear -ref %s',...
            sTimg,...
            fullfile(aap.directory_conventions.fsldir,'etc/flirtsch/ident.mat'),...
            spm_file(sTimg,'path',localroot,'prefix','r'),...
            fullfile(localroot,'tmp.nii')...
            ));
        
        % Rigid-body aligmnent to MNI with same image geometry
        aas_runfslcommand(aap,sprintf('flirt -in %s -ref %s -out %s -omat %s -dof 6',...
            Simg,...
            spm_file(sTimg,'path',localroot,'prefix','r'),...
            spm_file(Simg,'prefix','r'),...
            spm_file(Simg,'prefix','r','ext','mat')...
            ));
        
        %% Output
        % Clean-up
        delete(fullfile(localroot,'tmp.nii'));
        delete(spm_file(sTimg,'path',localroot,'prefix','r'));
        
        aap = aas_desc_outputs(aap,subj,'structural',spm_file(Simg,'prefix','r'));
end
end