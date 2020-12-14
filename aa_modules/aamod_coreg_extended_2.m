% AA module - extended coregistration
% Coregistration of structural to T1 template
% 1) Reorient MeanEPI to structural(output 'structural' from aamod_coreg_extended_1)
% 2) Coregister MeanEPI to structural (using xfm mat called
% 't1totemplate_xfm from aamod_coreg_extended_1)
% 3) Apply all (1 and 2) to epi.

function [aap,resp]=aamod_coreg_extended_2(aap,task,varargin)

resp='';

switch task
    case 'report' % [TA]
        domain = aap.tasklist.currenttask.domain;
        localpath = aas_getpath_bydomain(aap,domain,cell2mat(varargin));
        
        % Process streams
        [diagstream, mainstream] = process_streams(aap);
        % find out what streams we should normalise
        if isfield(aap.tasklist.currenttask.settings,'diagnostic') && isstruct(aap.tasklist.currenttask.settings.diagnostic)
            inds = aap.tasklist.currenttask.settings.diagnostic.streamind;
        else
            inds = 1:length(mainstream);
        end
        d = dir(fullfile(localpath,'diagnostic_aas_checkreg_*'));
        if isempty(d)
            if numel(varargin) > 1, diagstream = fullfile(aas_getsesspath(aap,varargin{1}, varargin{2}),'sessref.nii'); end
            aas_checkreg(aap,domain,cell2mat(varargin),diagstream,'structural');
            for m = inds
                if ~aas_stream_has_contents(aap,domain,cell2mat(varargin),mainstream{m}), continue; end
                aas_checkreg(aap,domain,cell2mat(varargin),mainstream{m},'structural');
            end
        end
        subj = varargin{1};
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
        flags = aap.spm.defaults.coreg;
        if isfield(aap.tasklist.currenttask.settings,'eoptions')
            fields = fieldnames(aap.tasklist.currenttask.settings.eoptions);
            for f = 1:numel(fields)
                if ~isempty(aap.tasklist.currenttask.settings.eoptions.(fields{f}))
                    flags.estimate.(fields{f}) = aap.tasklist.currenttask.settings.eoptions.(fields{f});
                end
            end
        end
        
        %% 0) Check that the templates and images we need exist!
        subj = varargin{1};
        
        % Process streams
        domain = aap.tasklist.currenttask.domain;
        [diagstream, mainstream, wbstream] = process_streams(aap);
        
        % Check local structural
        Simg = aas_getfiles_bystream(aap,subj,'structural');
        if size(Simg,1) > 1
            aas_log(aap, false, 'WARNING: Found more than 1 structural images, using first.');
            Simg = [deblank(Simg(1,:)),',1'];
        end
        
        % Look for mean functional
        mEPIimg = aas_getfiles_bystream_multilevel(aap,domain,cell2mat(varargin),diagstream);
        if numel(spm_vol(mEPIimg)) > 1
            aas_log(aap, false, 'WARNING: Found more than 1 mean functional images, using first.');
            mEPIimg = [deblank(mEPIimg(1,:)),',1'];
        end
        if numel(varargin) > 1
            V = spm_vol(mEPIimg); Y = spm_read_vols(V);
            V.fname = strtok(spm_file(mEPIimg,'path',aas_getpath_bydomain(aap,aap.tasklist.currenttask.domain,cell2mat(varargin)),'basename','sessref'),',');
            spm_write_vol(V,Y);
            mEPIimg = V.fname;
        end
        
        % Check local wholebrain EPI
        WBimg = '';
        if ~isempty(wbstream) && aas_stream_has_contents(aap,domain,cell2mat(varargin),wbstream)
            WBimg = aas_getfiles_bystream_multilevel(aap,domain,cell2mat(varargin),wbstream);
            if size(WBimg,1) > 1
                aas_log(aap, false, 'WARNING: Found more than 1 wholebrain images, using first.');
                WBimg = [deblank(WBimg(1,:)),',1'];
            end
            if numel(varargin) > 1
                copyfile(WBimg,spm_file(WBimg,'path',aas_getpath_bydomain(aap,aap.tasklist.currenttask.domain,cell2mat(varargin)),'basename','wholebrain'));
                WBimg = spm_file(WBimg,'path',aas_getpath_bydomain(aap,aap.tasklist.currenttask.domain,cell2mat(varargin)),'basename','wholebrain');
            end
        end
        
        % Look for xfm t1totemplate
        load(aas_getfiles_bystream(aap,subj,'t1totemplate_xfm'));
        
        aas_log(aap,false,sprintf(['\tto template realignment parameters:\n' ...
            '\tx: %0.4f   y: %0.4f   z: %0.4f   p: %0.4f   r: %0.4f   j: %0.4f'], ...
            xfm(1), xfm(2), xfm(3), xfm(4), xfm(5), xfm(6)))
        
        %% 1) Mean Functional to T1 template (reorient)
        % Set the new space for the mean functional
        for img = {mEPIimg WBimg}
            if isempty(img{1}), continue; end
            spm_get_space(img{1}, spm_matrix(xfm)\spm_get_space(img{1}));
        end

        %% 2) Mean Functional to Structural (coregister)

        if ~isempty(WBimg) % two-step coreg
            % Coregister wholebrain EPI to structural
            x1 = spm_coreg(spm_vol(Simg), ...
                spm_vol(WBimg), ...
                flags.estimate);
            % Set the new space for the wholebrain EPI
            spm_get_space(WBimg, spm_matrix(x1)\spm_get_space(WBimg));
            
            % Coregister mean EPI to wholebrain EPI
            flags.estimate.cost_fun = 'ncc'; % whithin-modality
            x2 = spm_coreg(spm_vol(deblank(WBimg)),spm_vol(mEPIimg),flags.estimate);
            % Set the new space for the mean EPI
            spm_get_space(mEPIimg, spm_matrix(x2)\spm_get_space(mEPIimg));
            x = x1 + x2;
        else
            % Coregister mean EPI to structural
            x = spm_coreg(spm_vol(Simg), ...
                spm_vol(mEPIimg), ...
                flags.estimate);
            % Set the new space for the mean EPI
            spm_get_space(mEPIimg, spm_matrix(x)\spm_get_space(mEPIimg));
        end
                
        aas_log(aap,false,sprintf(['\tmean EPI to structural realignment parameters:\n' ...
            '\tx: %0.4f   y: %0.4f   z: %0.4f   p: %0.4f   r: %0.4f   j: %0.4f'], ...
            x(1), x(2), x(3), x(4), x(5), x(6)))
        
        %% 3) Now apply this transformation to all the EPI images
        % The mean EPI will already be in the space required for the
        % individual EPIs. Hence, we can...
        
        % Again, get space of mean functional
        MM = spm_get_space(mEPIimg);
        
        % Locate all the EPIs we want to coregister
        for m = 1:numel(mainstream)
            if ~aas_stream_has_contents(aap,domain,cell2mat(varargin),mainstream{m}), continue; end
            EPIimg{m} = aas_getfiles_bystream(aap,domain,cell2mat(varargin),mainstream{m});
            excl = [];
            for e = 1:size(EPIimg{m},1)
                if ~strcmp(spm_file(EPIimg{m}(e,:),'ext'),'nii')
                    aas_log(aap,false,sprintf('WARNING: file %s is not a NIfTI --> skipping',EPIimg{m}(e,:)));
                    excl(end+1) = e;
                    continue; 
                end
                % Apply the space of the coregistered mean EPI to the
                % remaining EPIs (safest solution!)
                spm_get_space(deblank(EPIimg{m}(e,:)), MM);
            end
            EPIimg{m}(excl,:) = [];
        end
        
        %% Describe the outputs and Diagnostics
        
        if strcmp(aap.options.wheretoprocess,'localsingle')
            aas_checkreg(aap,domain,cell2mat(varargin),mEPIimg,'structural');
            for m = 1:numel(mainstream)
                aas_checkreg(aap,domain,cell2mat(varargin),mainstream{m},'structural');
            end
        end
        
        for m = 1:numel(mainstream)
            if ~aas_stream_has_contents(aap,domain,cell2mat(varargin),mainstream{m}), continue; end
            aap = aas_desc_outputs(aap,domain,cell2mat(varargin),mainstream{m},EPIimg{m});
        end
        if any(strcmp(aas_getstreams(aap,'output'),diagstream))
            aap = aas_desc_outputs(aap,domain,cell2mat(varargin),diagstream,mEPIimg);
        end
        
    case 'checkrequirements'
        in = aas_getstreams(aap,'input'); in(1:4) = []; % not for reference
        [stagename, index] = strtok_ptrn(aap.tasklist.currenttask.name,'_0');
        stageindex = sscanf(index,'_%05d');
        out = aap.tasksettings.(stagename)(stageindex).outputstreams.stream; if ~iscell(out), out = {out}; end
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

function [diagstream, mainstream, wbstream] = process_streams(aap)
[inpstreams, inpattr] = aas_getstreams(aap,'input');
outpstreams = aas_getstreams(aap,'output');

% select non-structural diagnostic stream
diagind = cellfun(@(x) isfield(x,'diagnostic') && x.diagnostic, inpattr); diagind(cell_index(inpstreams,'structural')) = false;
if any(diagind)
    diagstream = inpstreams{diagind};
else  % auto failed --> manual  
    if cell_index(inpstreams,'meanepi') % fMRI
        diagstream = 'meanepi';
    end
    if cell_index(inpstreams,'MTI_baseline') % MTI
        diagstream = 'MTI_baseline';
    end
end
diagstream = textscan(diagstream,'%s','delimiter','.'); diagstream = diagstream{1}{end};
diagind(cell_index(inpstreams,diagstream)) = true;

if cell_index(inpstreams,'wholebrain') % partial volume acquisition
    wbstream = inpstreams{cell_index(inpstreams,'wholebrain')};
    wbstream = textscan(wbstream,'%s','delimiter','.'); wbstream = wbstream{1}{end};
else
    wbstream = '';
end

% main: = out + ~diag
mainind = ~diagind;
[junk,ia] = intersect(inpstreams,outpstreams);
outind = false(size(mainind)); outind(ia) = true;
mainind = mainind & outind;
if any(mainind)
    mainstream = inpstreams(mainind);
else
    mainstream = {};
end
end