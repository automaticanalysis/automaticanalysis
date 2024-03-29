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
        [workstream, regstreams] = processStreams(aap,cell2mat(varargin));
        
        for w = workstream
            streamfn = aas_getfiles_bystream(aap,domain,cell2mat(varargin),w{1},'output');
            streamfn = strtok_ptrn(basename(streamfn(1,:)),'-0');
            fn = ['diagnostic_aas_checkreg_slices_' streamfn '_1.jpg'];
            ws = strsplit(w{1},'.'); ws = ws{end};
            if ~exist(fullfile(localpath,fn),'file')
                aas_checkreg(aap,domain,cell2mat(varargin),ws,regstreams{end});
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
        end
    case 'doit'
        %% Preapre
        subj = varargin{1};
        domain = aap.tasklist.currenttask.domain;
        domaindir = aas_getpath_bydomain(aap,domain,cell2mat(varargin));
        [workstream, regstreams] = processStreams(aap,cell2mat(varargin));

        for w = workstream
            images = cellstr(aas_getfiles_bystream_multilevel(aap, domain, cell2mat(varargin), w{1}));
            % make a working copy in the working folder to avoid interference
            tmp_images = spm_file(images,'path',domaindir,'suffix','_tmp');
            arrayfun(@(i) copyfile(images{i}, tmp_images{i}), 1:numel(images));
            images = tmp_images;
            
            %% Denorm
            if aas_stream_has_contents(aap,'inverse_deformation_field') % unified segment+normalise
                flags = aap.spm.defaults.normalise.write;
                [flags.bb, flags.vox] = spm_get_bbox(spm_vol(aas_getfiles_bystream_multilevel(aap,domain,cell2mat(varargin),regstreams{1})));
                flags.interp = aap.tasklist.currenttask.settings.interp;
                
                job.subj.def = cellstr(aas_getfiles_bystream(aap, subj, 'inverse_deformation_field'));
                job.subj.resample = images;
                job.woptions = flags;
                out = spm_run_norm(job);
            elseif aas_stream_has_contents(aap,'dartel_templatetomni_xfm') && aas_stream_has_contents(aap,'dartel_flowfield') % DARTEL
                ff = aas_getfiles_bystream(aap, subj, 'dartel_flowfield');
                if ~exist(fullfile(aas_getpath_bydomain(aap,domain,cell2mat(varargin)),[basename(ff) '.nii']),'file')
                    copyfile(ff,fullfile(aas_getpath_bydomain(aap,domain,cell2mat(varargin)),[basename(ff) '.nii']));
                end
                ff = fullfile(aas_getpath_bydomain(aap,domain,cell2mat(varargin)),basename([ff '.nii']));
                
                job.flowfields = {ff};
                job.images = images;
                job.interp = aap.tasklist.currenttask.settings.interp;
                job.K = 6;
                out = spm_dartel_invnorm(job);                
            elseif aas_stream_has_contents(aap,'freesurfer') % freesurfer
                fstemplate = fullfile(aas_getstudypath(aap),'fsaverage');
                if ~exist(fstemplate,'dir'), aas_shell(sprintf('ln -s %s/subjects/fsaverage %s',aap.directory_conventions.freesurferdir,fstemplate)); end
                setenv('SUBJECTS_DIR', aas_getstudypath(aap));
                for hemi = {'lh' 'rh'}
                    annot = images{cellfun(@(x) ~isempty(regexp(spm_file(x,'basename'),['^' hemi{1} '-.*'],'once')), images)};
                    FScommand = sprintf('mri_surf2surf --hemi %s --srcsubject fsaverage --trgsubject %s --sval-annot %s --tval %s',hemi{1},aas_getsubjname(aap,subj),annot,...
                        fullfile(aas_getsubjpath(aap,subj),'label',strrep(spm_file(annot,'filename'),[hemi{1} '-0'],hemi{1})));
                    aas_runFScommand(aap,FScommand);
                end
                out = [];
                out.files = cellstr(spm_file(annot,...
                    'path',aas_getsubjpath(aap,subj),...
                    'basename',strrep(spm_file(annot,'basename'),'rh-0.',''),...
                    'ext','nii'));
                FScommand = sprintf('mri_aparc2aseg --s %s --o %s --annot %s',aas_getsubjname(aap,subj),out.files{1},spm_file(out.files{1},'basename'));
                aas_runFScommand(aap,FScommand);
            else
                aas_log(aap,true,'No transformation secified')
            end
            
            %% (Decoreg +) Reslice
            % flags
            flags = aap.spm.defaults.coreg.write;
            flags.interp = aap.tasklist.currenttask.settings.interp;
            flags.which = [1 0];
            
            % Reslice
            if ~aas_getsetting(aap,'decoregextended')
                regstreams = regstreams(1);
            end
            imreslice = {aas_getfiles_bystream_multilevel(aap,domain,cell2mat(varargin),regstreams{end})};
            for r = 1:numel(out.files)
                imreslice{r+1} = spm_file(out.files{r},'path',aas_getpath_bydomain(aap,domain,cell2mat(varargin)));
                rwimgs{r} = spm_file(imreslice{r+1},'prefix',flags.prefix);
                if ~strcmp(out.files{r},imreslice{r+1}), movefile(out.files{r},imreslice{r+1}); end
                if numel(regstreams) == 2 % Decoreg + Reslice to final target
                    xfm = spm_get_space(aas_getfiles_bystream_multilevel(aap,domain,cell2mat(varargin),regstreams{2}))/...
                        spm_get_space(aas_getfiles_bystream_multilevel(aap,domain,cell2mat(varargin),regstreams{1}));
                    spm_get_space(imreslice{r+1},xfm*spm_get_space(imreslice{r+1}));
                end
            end
            spm_reslice(imreslice,flags);
            aap=aas_desc_outputs(aap,domain,cell2mat(varargin),w{1},rwimgs);
            
            %% Diag
            if strcmp(aap.options.wheretoprocess,'localsingle')
                aas_checkreg(aap,domain,cell2mat(varargin),w{1},regstreams{end});
            end
            
            %% Cleanup
            cellfun(@(f) delete(f), tmp_images);
        end
    case 'checkrequirements'
        in = processStreams(aap,cell2mat(varargin));
        [stagename, index] = strtok_ptrn(aap.tasklist.currenttask.name,'_0');
        stageindex = sscanf(index,'_%05d');
        out = aap.tasksettings.(stagename)(stageindex).outputstreams.stream; if ~iscell(out), out = {out}; end
        instream = cellfun(@(s) strsplit(s,'.'), in, 'UniformOutput', false); instream = cellfun(@(c) c(end), instream);
        for i = 1:numel(instream)
            if i <= numel(out)
                if ~strcmp(out{i},instream{i})
                    aap = aas_renamestream(aap,aap.tasklist.currenttask.name,out{i},instream{i},'output');
                    aas_log(aap,false,['INFO: ' aap.tasklist.currenttask.name ' output stream: ''' instream{i} '''']);
                end
            else
                aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'append',instream{i},'output');
                aas_log(aap,false,['INFO: ' aap.tasklist.currenttask.name ' output stream: ''' instream{i} '''']);
            end
        end
end
end

function [workstream, regstreams] = processStreams(aap,indices)
[inpstreams, streamattr] = aas_getstreams(aap,'input');
workstreamind = cellfun(@(a) ~isempty(a) && isfield(a,'isrenameable') && a.isrenameable, streamattr);
workstream = inpstreams(workstreamind); inpstreams(workstreamind) = [];

regstreamsind = cellfun(@(a) ~isempty(a) && isfield(a,'diagnostic') && a.diagnostic, streamattr);
regstreams = inpstreams(regstreamsind); toRemove = [];
for i = 1:numel(regstreams)
    if ~aas_stream_has_contents(aap,aap.tasklist.currenttask.domain,indices,regstreams{i}), toRemove(end+1) = i; end
end
regstreams(toRemove) = [];
end