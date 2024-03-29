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
                [p, f] = fileparts(imgpath); avipath = fullfile(p,[strrep(f(1:end-2),'slices','avi') '.avi']);
                if exist(avipath,'file'), aap=aas_report_addimage(aap,subj,avipath); end
                aap = aas_report_add(aap,subj,'</td></tr></table>');
            end
            % Study summary
            if all(cell2mat(varargin) == 1) % first
                stagename = regexp(aas_getstagetag(aap,aap.tasklist.currenttask.modulenumber),'_(?=0)','split'); stagename = stagename{1};
                stagename = [stagename aap.tasklist.currenttask.extraparameters.aap.directory_conventions.analysisid_suffix];
                aap = aas_report_add(aap,'reg',sprintf('<h2>Stage: %s</h2>',stagename));
            end
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
         [~,nam,ext] = fileparts(template);
        if ~exist(fullfile(localpath,[nam ext]),'file')
            copyfile(template,fullfile(localpath,[nam ext]));
        end
        template = fullfile(localpath,[nam ext]);
        
        % affine xfm
        streams=aas_getstreams(aap,'input');
        xfmi = cell_index(streams,'dartel_templatetomni_xfm');
        if xfmi && aas_stream_has_contents(aap,subj,streams{xfmi})
            
            verbose = aap.options.verbose;
            aap.options.verbose = -1;
            xfmfile = aas_getfiles_bystream_multilevel(aap,aap.tasklist.currenttask.domain,cell2mat(varargin),'dartel_templatetomni_xfm');
            aap.options.verbose = verbose;
            
            % Must be in a session folder
            if isempty(xfmfile)
                xfmfile = cellstr(spm_select('FPListRec',aas_getpath_bydomain(aap,aap.tasklist.currenttask.domain,cell2mat(varargin)),'^dartel_templatetomni_xfm.mat$'));
                xfmfile = xfmfile{1};
            end
            
            xfm = load(xfmfile); xfm = xfm.xfm;
            MMt = spm_get_space(template);
            mni.code = 'MNI152';
            mni.affine = xfm*MMt;
            [~,nam,~] = fileparts(template);
            save(fullfile(localpath,[nam '_2mni.mat']),'mni');            
        end
        
        % flow fields..
        job.data.subj.flowfield{1} = aas_getfiles_bystream(aap, subj, 'dartel_flowfield');
        job.template{1} = template;
        job.bb = nan(2,3);
        boundingBox = aas_getsetting(aap,'bb');
        if ~isempty(boundingBox)
            job.bb = reshape(boundingBox,2,3);
        end
        job.vox = aap.tasklist.currenttask.settings.vox;    % voxel size
        job.fwhm = aap.tasklist.currenttask.settings.fwhm;              % smoothing
        job.preserve = aap.tasklist.currenttask.settings.preserve;      % modulation
        
		% find out what streams we should normalise
        streams=aas_getstreams(aap,'output');
        for streamind=1:length(streams)
            if isstruct(streams{streamind}), streams{streamind} = streams{streamind}.CONTENT; end
            if ~aas_stream_has_contents(aap,streams{streamind}), continue; end            
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
            % delete any previous results because otherwise nifti write routine
            % doesn't save disc space when you reslice to a coarser voxel
            for c=1:size(P,1)
                [pth, fle, ext]=fileparts(P(c,:));
                previous = fullfile(pth,[aap.spm.defaults.normalise.write.prefix fle ext]);
                if exist(previous,'file') > 0
                    [s, w] = aas_shell(['rm ' previous],true); % quietly
                end
            end
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
            if ~isempty(aas_getsetting(aap,'PVE'))
                pve = aas_getsetting(aap,'PVE',streamind-1); % account for the skipped dartel_templatetomni_xfm
                if pve
                    for e = 1:size(wimgs,1)
                        V = spm_vol(deblank(wimgs(e,:)));
                        Y = spm_read_vols(V);
                        Y = Y >= pve;
                        nifti_write(deblank(wimgs(e,:)),Y,'Binarized',V)
                    end
                end
            end
            
            % describe outputs with diagnostoc
            aap=aas_desc_outputs(aap,index{:},streams{streamind},wimgs);
            if ~isfield(aap.tasklist.currenttask.settings,'diagnostic') ||...
                    (~isstruct(aap.tasklist.currenttask.settings.diagnostic) && aap.tasklist.currenttask.settings.diagnostic) ||...
                    (isstruct(aap.tasklist.currenttask.settings.diagnostic) && isfield(aap.tasklist.currenttask.settings.diagnostic,'streamind') && ((streamind - any(strcmp(streams,'dartel_templatetomni_xfm'))) == aap.tasklist.currenttask.settings.diagnostic.streamind))
                [inp, inpattr] = aas_getstreams(aap,'input');
                streamStruct = inp{cellfun(@(a) isfield(a,'diagnostic') && a.diagnostic, inpattr)};
                aas_checkreg(aap,aap.tasklist.currenttask.domain,cell2mat(varargin),streams{streamind},streamStruct);
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
        [in, inAttr] = aas_getstreams(aap,'input'); 
        [in, indIn] = setdiff(in,{'dartel_template' 'dartel_templatetomni_xfm' 'dartel_flowfield'},'stable'); inAttr = inAttr(indIn);
        % do not consider structural if it is a diagnostic stream
        if any(contains(in,'structural')) && isfield(inAttr{contains(in,'structural')},'diagnostic') && inAttr{contains(in,'structural')}.diagnostic
            in(contains(in,'structural')) = [];
        end
        out = aas_getstreams(aap,'output'); out = setdiff(out,{'dartel_templatetomni_xfm'},'stable');
        for s = 1:numel(in)
            instream = textscan(in{s},'%s','delimiter','.'); instream = instream{1}{end};
            if s <= numel(out)
                if ~strcmp(out{s},instream)
                    aap = aas_renamestream(aap,aap.tasklist.currenttask.name,out{s},instream,'output');
                    aas_log(aap,false,['INFO: ' aap.tasklist.currenttask.name ' output stream: ''' instream '''']);
                end
            else
                aap = aas_renamestream(aap,aap.tasklist.currenttask.name,'append',instream,'output');
                aas_log(aap,false,['INFO: ' aap.tasklist.currenttask.name ' output stream: ''' instream '''']);
            end
        end
        
end
end
