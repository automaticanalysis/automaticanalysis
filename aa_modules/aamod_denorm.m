function [aap, resp]=aamod_denorm(aap, task, varargin)
% Inputstream order
%   1. coreg_exteded.main - subject domain (optional)
%   2. native.session - session domain
%   3. inverse_deformation_field
%   4. work stream in MNI space

resp='';

% possible tasks 'doit','report','checkrequirements'
switch task
    case 'report'
        subj = varargin{1};
        domain = aap.tasklist.currenttask.domain;
        localpath = aas_getpath_bydomain(aap,domain,cell2mat(varargin));
        inpstreams = aas_getstreams(aap,'input');
        workstream = inpstreams{end};
        images = cellstr(aas_getfiles_bystream_multilevel(aap, domain, cell2mat(varargin), workstream));
        regstreams = inpstreams(1:end-3); toRemove = [];
        for i = 1:numel(regstreams)
            if ~aas_stream_has_contents(aap,domain,cell2mat(varargin),regstreams{i}), toRemove(end+1) = i; end
        end
        regstreams(toRemove) = [];

        streamfn = aas_getfiles_bystream(aap,domain,cell2mat(varargin),workstream,'output');
        streamfn = strtok_ptrn(basename(streamfn(1,:)),'-0');
        fn = ['diagnostic_aas_checkreg_slices_' streamfn '_1.jpg'];
        workstream = strsplit(workstream,'.'); workstream = workstream{end};
        if ~exist(fullfile(localpath,fn),'file')
            aas_checkreg(aap,domain,cell2mat(varargin),workstream,regstreams{end});
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
        subj = varargin{1};
        domain = aap.tasklist.currenttask.domain;
        inpstreams = aas_getstreams(aap,'input');
        workstream = inpstreams{end};
        images = cellstr(aas_getfiles_bystream_multilevel(aap, domain, cell2mat(varargin), workstream));
        regstreams = inpstreams(1:end-2); toRemove = [];
        for i = 1:numel(regstreams)
            if ~aas_stream_has_contents(aap,domain,cell2mat(varargin),regstreams{i}), toRemove(end+1) = i; end
        end
        regstreams(toRemove) = [];

        %% Denorm
        flags = aap.spm.defaults.normalise.write;
        [flags.bb, flags.vox] = spm_get_bbox(spm_vol(aas_getfiles_bystream_multilevel(aap,domain,cell2mat(varargin),regstreams{1})));
        flags.interp = aap.tasklist.currenttask.settings.interp;
        
        job.subj.def = cellstr(aas_getfiles_bystream(aap, subj, 'inverse_deformation_field'));
        job.subj.resample = images;
        job.woptions = flags;
        out = spm_run_norm(job);
        
        %% (Decoreg +) Reslice
        % flags
        flags = aap.spm.defaults.coreg.write;
        flags.interp = aap.tasklist.currenttask.settings.interp;
        flags.which = [1 0];

        % Reslice
        imreslice = {aas_getfiles_bystream_multilevel(aap,domain,cell2mat(varargin),regstreams{end})};
        for r = 1:numel(out.files)
            imreslice{r+1} = spm_file(out.files{r},'path',aas_getpath_bydomain(aap,domain,cell2mat(varargin)));
            rwimgs{r} = spm_file(imreslice{r+1},'prefix',flags.prefix);
            movefile(out.files{r},imreslice{r+1});
            if numel(regstreams) == 2 % Decoreg + Reslice to final target
                xfm = spm_get_space(aas_getfiles_bystream_multilevel(aap,domain,cell2mat(varargin),regstreams{2}))/...
                    spm_get_space(aas_getfiles_bystream_multilevel(aap,domain,cell2mat(varargin),regstreams{1}));
                spm_get_space(imreslice{r+1},xfm*spm_get_space(imreslice{r+1}));
            end
        end
        spm_reslice(imreslice,flags);
        aap=aas_desc_outputs(aap,domain,cell2mat(varargin),workstream,rwimgs);
        
        %% Diag
        if strcmp(aap.options.wheretoprocess,'localsingle')
            aas_checkreg(aap,domain,cell2mat(varargin),workstream,regstreams{end});
        end
    case 'checkrequirements'
        in =  aas_getstreams(aap,'input'); in = in{end}; % last stream
        [stagename, index] = strtok_ptrn(aap.tasklist.currenttask.name,'_0');
        stageindex = sscanf(index,'_%05d');
        out = aap.tasksettings.(stagename)(stageindex).outputstreams.stream; if iscell(out), out = out{1}; end
        instream = textscan(in,'%s','delimiter','.'); instream = instream{1}{end};
        if ~strcmp(out,instream)
            aap = aas_renamestream(aap,aap.tasklist.currenttask.name,out,instream,'output');
            aas_log(aap,false,['INFO: ' aap.tasklist.currenttask.name ' output stream: ''' instream '''']);
        end
end