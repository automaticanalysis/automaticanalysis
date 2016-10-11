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
		streams=aas_getstreams(aap,'output');
        if isfield(aap.tasklist.currenttask.settings,'diagnostic') && isstruct(aap.tasklist.currenttask.settings.diagnostic)
            inds = aap.tasklist.currenttask.settings.diagnostic.streamind;
        else
            inds = 1:length(streams);
        end
        % determine normalised struct
        struct = aas_getfiles_bystream_dep(aap,'subject',varargin{1},'structural');
        sname = basename(struct);
        struct = struct((sname(:,1)=='w') | (sname(:,2)=='w'),:);
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
            aap=aas_report_addimage(aap,'reg',fullfile(localpath,fdiag(1).name));
        end
    case 'doit'
        subj = varargin{1};
        if nargin == 4, sess = varargin{2}; end
        
        % Is session specified in task header?
        if (isfield(aap.tasklist.currenttask.settings,'session'))
            sess = aap.tasklist.currenttask.settings.session;
        end        
        
        if aas_stream_has_contents(aap,'normalisation_seg_sn')
            try flags = aap.spm.defaults.old.normalise.write; catch, flags = aap.spm.defaults.normalise.write; end
            % get sn mat file from normalisation
            trans = aas_getfiles_bystream(aap,subj,'normalisation_seg_sn');
        elseif aas_stream_has_contents(aap,'forward_deformation_field')
            flags = aap.spm.defaults.normalise.write;
            trans = aas_getfiles_bystream(aap,subj,'forward_deformation_field');
        else
            aas_log(aap,true,'ERROR: No transformation specified!')
        end
        
        flags.vox = aas_getsetting(aap,'vox');
        boundingBox = aas_getsetting(aap,'bb');
        if ~isempty(boundingBox), flags.bb = reshape(boundingBox,2,3); end
        interp = aas_getsetting(aap,'interp');
        if ~isempty(interp), flags.interp = interp; end
        
		% find out what streams we should normalise
        streams=aap.tasklist.currenttask.outputstreams.stream;
        for streamind=1:length(streams)
            % Image to reslice
            if isstruct(streams{streamind}), streams{streamind} = streams{streamind}.CONTENT; end
            P = aas_getfiles_bystream(aap,aap.tasklist.currenttask.domain,cell2mat(varargin),streams{streamind});
            
            % Ignore .hdr files from this list...
            ishdr = cell_index(cellstr(P),'.hdr');
            if any(ishdr), P(ishdr,:) = []; end
            
            % delete previous because otherwise nifti write routine doesn't
            % save disc space when you reslice to a coarser voxel
            for c=1:size(P,1)
                aas_shell(['rm ' spm_file(P(c,:),'prefix',flags.prefix)],true); % quietly
            end
            
            if aas_stream_has_contents(aap,'normalisation_seg_sn')
                spm_write_sn(P,trans,flags);
            elseif aas_stream_has_contents(aap,'forward_deformation_field')
                switch spm('ver')
                    case 'SPM8'
                        job.ofname = '';
                        job.fnames = cellstr(P);
                        job.savedir.saveusr{1} = aas_getpath_bydomain(aap,aap.tasklist.currenttask.domain,cell2mat(varargin));
                        job.interp = flags.interp;
                        job.comp{1}.def = cellstr(trans);
                        spm_defs(job);
                    case {'SPM12b' 'SPM12'}
                        job.subj.def = cellstr(trans);
                        job.subj.resample = cellstr(P);
                        job.woptions = flags; 
                        spm_run_norm(job);
                    otherwise
                        aas_log(aap, true, sprintf('%s requires SPM8 or later.', mfilename));
                end
            end
            
            % outputs
            wimgs=[];
            for fileind=1:size(P,1)
                % overwrite input with output if specified (e.g. for contrasts)
                if isfield(aap.tasklist.currenttask.settings.outputstreams,'preservefilename') && ...
                        aap.tasklist.currenttask.settings.outputstreams.preservefilename
                    movefile(spm_file(P(fileind,:),'prefix',flags.prefix),P(fileind,:));
                    wimgs = strvcat(wimgs,P(fileind,:));
                else
                    wimgs = strvcat(wimgs,spm_file(P(fileind,:),'prefix',flags.prefix));
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
            aap=aas_desc_outputs(aap,aap.tasklist.currenttask.domain,cell2mat(varargin),streams{streamind},wimgs);
            if strcmp(aap.options.wheretoprocess,'localsingle') && any(strcmp(streams, 'structural'))
                aas_checkreg(aap,aap.tasklist.currenttask.domain,cell2mat(varargin),streams{streamind},struct);
            end
        end
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
end