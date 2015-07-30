% AA module - write normalised input useing DARTEL
% [aap,resp]=aamod_norm_write_dartel(aap,task,subj,sess)
% Rhodri Cusack MRC CBU Cambridge Jan 2006-Aug 2007
% Tibor Auer MRC CBU Cambridge 2012-2013

function [aap,resp]=aamod_norm_write_dartel(aap,task,varargin)

resp='';

switch task
    case 'report' % [TA]
        subj = varargin{1};
        if nargin == 4
            sess = varargin{2};
            localpath = aas_getpath_bydomain(aap,aap.tasklist.currenttask.domain,[subj,sess]);
        else % subject
            localpath = aas_getpath_bydomain(aap,'subject',subj);
        end
        
        % find out what streams we should normalise
		streams=aas_getstreams(aap,'output');
        if isfield(aap.tasklist.currenttask.settings,'diagnostic') && isstruct(aap.tasklist.currenttask.settings.diagnostic)
            inds = aap.tasklist.currenttask.settings.diagnostic.streamind;
        else
            inds = 1:length(streams);
        end
        for streamind = inds
            if strcmp(streams{streamind}, 'dartel_templatetomni_xfm'), continue, end
            streamfn = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,cell2mat(varargin),streams{streamind},'output');
            streamfn = streamfn(1,:);
            streamfn = strtok_ptrn(basename(streamfn),'-0');
            fn = ['diagnostic_aas_checkreg_slices_' streamfn '_1.jpg'];
            if ~exist(fullfile(localpath,fn),'file')
                aas_checkreg(aap,aap.tasklist.currenttask.domain,cell2mat(varargin),streams{streamind},'structural');
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
            % Study summary
            aap = aas_report_add(aap,'reg',...
                ['Subject: ' basename(aas_getsubjpath(aap,subj)) '; Session: ' aas_getdirectory_bydomain(aap,aap.tasklist.currenttask.domain,varargin{end}) ]);
            aap=aas_report_addimage(aap,'reg',fullfile(localpath,fdiag(1).name));
        end
    case 'doit'
        index{1} = aap.tasklist.currenttask.domain;
        subj = varargin{1};
        index{2}(1) = subj;
        if nargin == 4
            sess = varargin{2}; 
            index{2}(2) = sess;
        end
        
        % Is session specified in task header?
        if (isfield(aap.tasklist.currenttask.settings,'session'))
            sess = aap.tasklist.currenttask.settings.session;
            index{2}(2) = sess;
        end
        
        localpath = aas_getpath_bydomain(aap,index{:});
        
        % set up job
        % template
        template = aas_getfiles_bystream(aap, 'dartel_template');
         [pth,nam,ext] = fileparts(template);
        if ~exist(fullfile(localpath,[nam ext]),'file')
            copyfile(template,fullfile(localpath,[nam ext]));
        end
        template = fullfile(localpath,[nam ext]);
        
        % affine xfm
        streams=aas_getstreams(aap,'input');
        xfmi = cell_index(streams,'dartel_templatetomni_xfm');
        if xfmi && aas_stream_has_contents(aap,subj,streams{xfmi})
            xfm = load(aas_getfiles_bystream(aap,subj,'dartel_templatetomni_xfm')); xfm = xfm.xfm;
            MMt = spm_get_space(template);
            mni.code = 'MNI152';
            mni.affine = xfm*MMt;
            [pth,nam,ext] = fileparts(template);
            save(fullfile(localpath,[nam '_2mni.mat']),'mni');            
        end
        
        % flow fields..
        job.data.subj.flowfield{1} = aas_getfiles_bystream(aap, subj, 'dartel_flowfield');
        job.template{1} = template;
        job.bb = nan(2,3);
        job.vox = aap.tasklist.currenttask.settings.vox;    % voxel size
        job.fwhm = aap.tasklist.currenttask.settings.fwhm;              % smoothing
        job.preserve = aap.tasklist.currenttask.settings.preserve;      % modulation
        
		% find out what streams we should normalise
        streams=aap.tasklist.currenttask.outputstreams.stream;        
        for streamind=1:length(streams)
            if ~aas_stream_has_contents(aap,streams{streamind}), continue; end
            if isstruct(streams{streamind}), streams{streamind} = streams{streamind}.CONTENT; end
            if strcmp(streams{streamind},'dartel_templatetomni_xfm'), continue; end % skip            
            imgs = [];
            % Image to reslice
            if exist('sess','var')
                P = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,[subj,sess],streams{streamind});
            else
                P = aas_getfiles_bystream(aap,subj,streams{streamind});
            end
            % exclude image already normalised
            f = basename(P);
            P = P(f(:,1) ~= 'w',:);
            imgs = strvcat(imgs, P);
            % delete previous because otherwise nifti write routine doesn't
            % save disc space when you reslice to a coarser voxel
            for c=1:size(P,1)
                [pth fle ext]=fileparts(P(c,:));
                [s w] = aas_shell(['rm ' fullfile(pth,[aap.spm.defaults.normalise.write.prefix fle ext])],true); % quietly
            end;
            job.data.subj.images = cellstr(imgs);
            
            aas_log(aap, false, sprintf('Running with %s...', which('spm_dartel_norm_fun')));
            spm_dartel_norm_fun(job);
            
            % describe outputs (differ depending on modulation)
            prefix = 'w';
            if aap.tasklist.currenttask.settings.preserve, prefix = ['m' prefix]; end
            if aap.tasklist.currenttask.settings.fwhm, prefix = ['s' prefix]; end
        
            wimgs=[];
            for ind=1:length(job.data.subj.images)
                [pth, nme, ext] = fileparts(job.data.subj.images{ind});
                % overwrite input with output if specified (e.g. for contrasts)
                if isfield(aap.tasklist.currenttask.settings.outputstreams,'preservefilename') && ...
                        aap.tasklist.currenttask.settings.outputstreams.preservefilename
                    movefile(fullfile(pth,[prefix nme ext]),job.data.subj.images{ind});
                    wimgs = strvcat(wimgs,job.data.subj.images{ind});
                else
                    wimgs = strvcat(wimgs,fullfile(pth,[prefix nme ext]));
                end
            end
            
            % binarise if specified
            if isfield(aap.tasklist.currenttask.settings,'PVE') && ~isempty(aap.tasklist.currenttask.settings.PVE)
                for e = 1:size(wimgs,1)
                    inf = spm_vol(deblank(wimgs(e,:)));
                    Y = spm_read_vols(inf);
                    Y = Y>=aap.tasklist.currenttask.settings.PVE;
                    nifti_write(deblank(wimgs(e,:)),Y,'Binarized',inf)
                end
            end
            
            % describe outputs with diagnostoc
            aap=aas_desc_outputs(aap,index{:},streams{streamind},wimgs);
            if strcmp(aap.options.wheretoprocess,'localsingle')
                aas_checkreg(aap,index{:},streams{streamind},'structural');
            end
            if ~exist('xfm','var')
                MMt = spm_get_space(template);
                MMm = load(fullfile(localpath,[nam '_2mni.mat']));
                MMm = MMm.mni.affine;
                xfm = MMm/MMt;
            end
            save(fullfile(localpath,'dartel_templatetomni_xfm'),'xfm')
            aap = aas_desc_outputs(aap, index{:}, 'dartel_templatetomni_xfm',...
                fullfile(localpath,'dartel_templatetomni_xfm.mat'));
        end
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
end