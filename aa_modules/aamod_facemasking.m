function [aap,resp]=aamod_facemasking(aap,task,varargin)
resp='';

switch task
    case 'report'
    case 'doit'
        % init
        indices = cell2mat(varargin);        
        streams=aas_getstreams(aap,'input');
        
        % obtain input filename(s)
        infname = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,indices,streams{1});
        
        % do stuff
        % Setup
        wdir = fullfile(aas_getpath_bydomain(aap,aap.tasklist.currenttask.domain,indices),aap.directory_conventions.structdirname,'maskface');
        aas_makedir(aap,wdir);
        
        % Reslice input image
        % - resize template
        tmpl = fullfile(aap.directory_conventions.spmdir, aap.directory_conventions.T1template);
        copyfile(tmpl,spm_file(tmpl,'path',wdir))
        [junk,res] = spm_get_bbox(infname);
        resize_img(spm_file(tmpl,'path',wdir), abs(res), [nan nan nan; nan nan nan])        
        % - reslice input
        flags = aap.spm.defaults.coreg;
        flags.write.which = [1 0]; % do not (re)write target and mean
        flags.write.interp = 1; % spm_reslice default
        spm_reslice({spm_file(tmpl,'path',wdir,'prefix','r') infname},flags.write)
        infname = spm_file(infname,'prefix',flags.write.prefix);
        
        FslOT = aap.directory_conventions.fsloutputtype;
        aap.directory_conventions.fsloutputtype = 'NIFTI_PAIR';
        ENV = {...
            'MASKFACE_HOME', aap.directory_conventions.FaceMaskingdir;...
            'PATH', sprintf('${MASKFACE_HOME}/bin:%s/bin:${PATH}',matlabroot) ...
            };
        cmd = sprintf('fslmaths %s %s;cd %s;mask_face %s -a',...
            infname,...
            fullfile(wdir,spm_file(infname,'basename')),...
            fileparts(wdir),...
            fullfile(wdir,spm_file(infname,'basename')));
        aas_runfslcommand(aap,cmd,ENV);
        
        % Output
        % - write out .nii
        outfname = spm_file(infname,'prefix','defaced_');
        aap.directory_conventions.fsloutputtype = FslOT;
        aas_runfslcommand(aap,sprintf('fslmaths %s %s',...
            fullfile(wdir,[spm_file(infname,'basename') '_full_normfilter.hdr']),...
            outfname));
        % - correct space information (FSL-SPM issue)
        spm_get_space(outfname,spm_get_space(infname));
        
        aap=aas_desc_outputs(aap,aap.tasklist.currenttask.domain,indices,['defaced_' streams{1}],outfname);
    case 'checkrequirements'
        if ~strcmp(aap.internal.inputstreamsources{aap.tasklist.currenttask.modulenumber}.stream(1).sourcestagename,'aamod_coreg_extended_1')
            aas_log(aap,false,'WARNING: FaceMasking requires image to be coregistered into MNI space (e.g. aamod_coreg_extended_1)')
        end
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end