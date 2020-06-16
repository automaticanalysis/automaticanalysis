function [aap, resp] = aamod_fsl_reorienttoMNI(aap,task,varargin)

resp = '';

switch task
    case 'doit'
        [junk, SPMtool] = aas_cache_get(aap,'spm');
        
        %% Init
        subj = varargin{1};
        switch aap.tasklist.currenttask.domain
            case 'subject' % structural
                localroot = fullfile(aas_getsubjpath(aap,subj),aap.directory_conventions.structdirname);
            case 'session' % fmri
                localroot = aas_getsesspath(aap,subj,varargin{2});
        end
        
        % Get the T1 template
        sTimg = aap.directory_conventions.T1template;
        if ~exist(sTimg,'file') % try in SPM
            if sTimg(1) ~= '/', sTimg = fullfile(SPMtool.toolPath,sTimg); end
        else
            sTimg = which(sTimg);
        end
        if ~exist(sTimg,'file')
            aas_log(aap, true, sprintf('Couldn''t find template T1 image %s.', sTimg));
        end  
        
        % Get structural
        Stream = aas_getstreams(aap,'input');
        Simg = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,cell2mat(varargin),Stream{1});        
        
        %% Run
        % Obtain image geometry
        junk = spm_vol(sTimg);
        dim = abs(junk.dim);
        [junk, vox] = spm_get_bbox(sTimg); vox = abs(vox);
        bbox = dim.*vox;
        
        junk = spm_vol(Simg); junk = junk(1);
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
        
        aap = aas_desc_outputs(aap,aap.tasklist.currenttask.domain,cell2mat(varargin),Stream{1},spm_file(Simg,'prefix','r'));
    case 'checkrequirements'
        if ~aas_cache_get(aap,'spm'), aas_log(aap,true,'SPM is not found'); end
end
end