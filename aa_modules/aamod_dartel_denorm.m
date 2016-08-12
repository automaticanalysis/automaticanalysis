function [aap, resp]=aamod_dartel_denorm(aap, task, subj, sess)
% Inputstream order
%   1. coreg_exteded.main - subject domain (optional)
%   2. native.session - session domain
%   3. dartel_templatetomni_xfm
%   4. flowfield
%   5. work stream in MNI space

resp='';

% possible tasks 'doit','report','checkrequirements'
switch task
    case 'report'
        domain = aap.tasklist.currenttask.domain;
        localpath = aas_getpath_bydomain(aap,domain,[subj,sess]);
        inpstreams = aas_getstreams(aap,'input');
        workstream = inpstreams{end};
        images = cellstr(aas_getfiles_bystream_multilevel(aap, domain, [subj sess], workstream));
        regstreams = inpstreams(1:end-3); toRemove = [];
        for i = 1:numel(regstreams)
            if ~aas_stream_has_contents(aap,domain,[subj sess],regstreams{i}), toRemove(end+1) = i; end
        end
        regstreams(toRemove) = [];

        streamfn = aas_getfiles_bystream(aap,domain,[subj sess],['native_' workstream],'output');
        streamfn = strtok_ptrn(basename(streamfn(1,:)),'-0');
        fn = ['diagnostic_aas_checkreg_slices_' streamfn '_1.jpg'];
        if ~exist(fullfile(localpath,fn),'file')
            aas_checkreg(aap,domain,[subj,sess],['native_' workstream],regstreams{end});
        end
        % Single-subjetc
        fdiag = dir(fullfile(localpath,'diagnostic_*.jpg'));
        for d = 1:numel(fdiag)
            aap = aas_report_add(aap,subj,'<table><tr><td>');
            imgpath = fullfile(localpath,fdiag(d).name);
            aap=aas_report_addimage(aap,subj,imgpath);
            [p f] = fileparts(imgpath); avipath = fullfile(p,[strrep(f(1:end-2),'slices','avi') '.avi']);
            if exist(avipath,'file'), aap=aas_report_addimage(aap,subj,avipath); end
            aap = aas_report_add(aap,subj,'</td></tr></table>');
        end
        
    case 'doit'
        %% Preapre
        domain = aap.tasklist.currenttask.domain;
        localpath = aas_getpath_bydomain(aap,domain,[subj,sess]);
        inpstreams = aas_getstreams(aap,'input');
        workstream = inpstreams{end};
        images = cellstr(aas_getfiles_bystream_multilevel(aap, domain, [subj sess], workstream));
        regstreams = inpstreams(1:end-3); toRemove = [];
        for i = 1:numel(regstreams)
            if ~aas_stream_has_contents(aap,domain,[subj sess],regstreams{i}), toRemove(end+1) = i; end
        end
        regstreams(toRemove) = [];

        % Affine Template2mni
        xfm = load(aas_getfiles_bystream_multilevel(aap, domain, [subj sess], 'dartel_templatetomni_xfm'));
        for i = 1:numel(images)
            if ~exist(fullfile(aas_getpath_bydomain(aap,domain,[subj,sess]),[basename(images{i}) '.nii']),'file')
                copyfile(images{i},fullfile(aas_getpath_bydomain(aap,domain,[subj,sess]),[basename(images{i}) '.nii']));
            end
            images{i} = fullfile(aas_getpath_bydomain(aap,domain,[subj,sess]),basename([images{i} '.nii']));
            spm_get_space(images{i},xfm.xfm\spm_get_space(images{i}));
        end
        
        % Flowfield
        ff = aas_getfiles_bystream(aap, subj, 'dartel_flowfield');
        if ~exist(fullfile(aas_getpath_bydomain(aap,domain,[subj,sess]),[basename(ff) '.nii']),'file')
            copyfile(ff,fullfile(aas_getpath_bydomain(aap,domain,[subj,sess]),[basename(ff) '.nii']));
        end
        ff = fullfile(aas_getpath_bydomain(aap,domain,[subj,sess]),basename([ff '.nii']));
        
        %% Denorm
        job.flowfields = {ff};
        job.images = images;
        job.interp = aap.tasklist.currenttask.settings.interp;
        job.K = 6;
        spm_dartel_invnorm(job);
        
        %% (Decoreg +) Reslice
        % flags
        defaults = aap.spm.defaults;
        flags = defaults.coreg.write;
        flags.interp = aap.tasklist.currenttask.settings.interp;
        flags.which = [1 0];

        % Reslice
        imreslice{1} = aas_getfiles_bystream_multilevel(aap,domain,[subj sess],regstreams{end});
        for r = 1:numel(job.images)
            imreslice{r+1} = fullfile(aas_getpath_bydomain(aap,domain,[subj sess]),['w' basename(job.images{r}) '.nii']);
            rwimgs{r} = fullfile(aas_getpath_bydomain(aap,domain,[subj sess]),[flags.prefix 'w' basename(job.images{r}) '.nii']);
            movefile(fullfile(fileparts(job.flowfields{1}),['w' basename(job.images{r}) '_' basename(job.flowfields{1}) '.nii']),imreslice{r+1});
            if numel(regstreams) == 2 % Decoreg + Reslice to final target
                xfm = spm_get_space(aas_getfiles_bystream_multilevel(aap,domain,[subj sess],regstreams{2}))/...
                    spm_get_space(aas_getfiles_bystream_multilevel(aap,domain,[subj sess],regstreams{1}));
                spm_get_space(imreslice{r+1},xfm*spm_get_space(imreslice{r+1}))
            end
        end
        spm_reslice(imreslice,flags);
        aap=aas_desc_outputs(aap,domain,[subj,sess],['native_' workstream],rwimgs);
        
        %% Diag
        if strcmp(aap.options.wheretoprocess,'localsingle')
            aas_checkreg(aap,domain,[subj,sess],['native_' workstream],regstreams{end});
        end
end