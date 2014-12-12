% AA module - write normalised EPIs
% [aap,resp]=aamod_norm_write(aap,task,subj,sess)
% Rhodri Cusack MRC CBU Cambridge Jan 2006-Aug 2007
% Resamples EPIs using *_seg_sn.mat file [if present] or *_sn.mat file
% Changed domain to once per session for improved performance when parallel
% Tibor Auer MRC CBU Cambridge 2012-2013

function [aap,resp]=aamod_norm_write(aap,task,varargin)

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
		streams=aas_getstreams(aap,'out');
        if isfield(aap.tasklist.currenttask.settings,'diagnostic') && isstruct(aap.tasklist.currenttask.settings.diagnostic)
            inds = aap.tasklist.currenttask.settings.diagnostic.streamind;
        else
            inds = 1:length(streams);
        end
        % determine normalised struct
        struct = aas_getfiles_bystream(aap,'subject',varargin{1},'structural');
        sname = basename(struct);
        struct = struct((sname(:,1)=='w'),:);
        for streamind = inds
            streamfn = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,cell2mat(varargin),streams{streamind},'output');
            streamfn = streamfn(1,:);
            streamfn = strtok_ptrn(basename(streamfn),'-0');
            fn = ['diagnostic_aas_checkreg_slices_' streamfn '_1.jpg'];
            if ~exist(fullfile(localpath,fn),'file')
                aas_checkreg(aap,aap.tasklist.currenttask.domain,cell2mat(varargin),streams{streamind},struct);
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
            aap=aas_report_addimage(aap,'reg',fullfile(localpath,fn));
        end
    case 'doit'
        subj = varargin{1};
        if nargin == 4, sess = varargin{2}; end
        
        % Is session specified in task header?
        if (isfield(aap.tasklist.currenttask.settings,'session'))
            sess = aap.tasklist.currenttask.settings.session;
        end        
        voxelSize = aap.tasklist.currenttask.settings.vox; % in case we want something other than default voxel size
        
        % get sn mat file from normalisation
        matname = aas_getfiles_bystream(aap,subj,'normalisation_seg_sn');
        
		% find out what streams we should normalise
        streams=aap.tasklist.currenttask.outputstreams.stream;
        
        for streamind=1:length(streams)
            imgs = [];
            % Image to reslice
            if isstruct(streams{streamind}), streams{streamind} = streams{streamind}.CONTENT; end
            if exist('sess','var')
                P = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,[subj,sess],streams{streamind});
            else
                P = aas_getfiles_bystream(aap,subj,streams{streamind});
            end
            
            % exclude image already normalised
            f = basename(P);
            P = P(f(:,1) ~= aap.spm.defaults.normalise.write.prefix,:);
            
            imgs = strvcat(imgs, P);
            
            % delete previous because otherwise nifti write routine doesn't
            % save disc space when you reslice to a coarser voxel
            for c=1:size(P,1)
                [pth fle ext]=fileparts(P(c,:));
                [s w] = aas_shell(['rm ' fullfile(pth,[aap.spm.defaults.normalise.write.prefix fle ext])],true); % quietly
            end;
            
            
            % set defaults
            flags = aap.spm.defaults.normalise.write;
            flags.vox = voxelSize;
            
            % now write normalised
            if ~isempty(imgs)
                % Ignore .hdr files from this list...
                imgsGood = imgs;
                for n = size(imgsGood,1):-1:1
                    if ~isempty(strfind(imgsGood(n,:), '.hdr'))
                        imgsGood(n,:) = [];
                    end
                end
                spm_write_sn(imgsGood,matname,aap.spm.defaults.normalise.write);
            end
            
            wimgs=[];
            
            % describe outputs
            for fileind=1:size(imgs,1)
                [pth, nme, ext] = fileparts(imgs(fileind,:));
                % overwrite input with output if specified (e.g. for contrasts)
                if isfield(aap.tasklist.currenttask.settings.outputstreams,'preservefilename') && ...
                        aap.tasklist.currenttask.settings.outputstreams.preservefilename
                    movefile(fullfile(pth,[aap.spm.defaults.normalise.write.prefix nme ext]),imgs(fileind,:));
                    wimgs = strvcat(wimgs,imgs(fileind,:));
                else
                    wimgs = strvcat(wimgs,fullfile(pth,[aap.spm.defaults.normalise.write.prefix nme ext]));
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
            
            % describe outputs with diagnostic
            % determine normalised struct
            struct = aas_getfiles_bystream(aap,'subject',varargin{1},'structural');
            sname = basename(struct);
            struct = struct((sname(:,1)=='w'),:);
            if (exist('sess','var'))
                aap=aas_desc_outputs(aap,aap.tasklist.currenttask.domain,[subj,sess],streams{streamind},wimgs);
                if strcmp(aap.options.wheretoprocess,'localsingle') && ismember(streams, 'structural')
                    aas_checkreg(aap,aap.tasklist.currenttask.domain,[subj,sess],streams{streamind},struct);
                end
            else
                aap=aas_desc_outputs(aap,subj,streams{streamind},wimgs);
                if strcmp(aap.options.wheretoprocess,'localsingle') && ismember(streams, 'structural')
                    aas_checkreg(aap,subj,streams{streamind},struct);
                end
            end
        end
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
end