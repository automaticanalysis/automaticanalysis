% AA module - extended coregistration
% Coregistration of structural to T1 template
% 1) Reorient MeanEPI to structural(output 'structural' from aamod_coreg_extended_1)
% 2) Coregister MeanEPI to structural (using xfm mat called
% 't1totemplate_xfm from aamod_coreg_extended_1)
% 3) Apply all (1 and 2) to epi.

function [aap,resp]=aamod_epireg_extended_2(aap,task,varargin)

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
            [p, f] = fileparts(imgpath); avipath = fullfile(p,[strrep(f(1:end-2),'slices','avi') '.avi']);
            if exist(avipath,'file'), aap=aas_report_addimage(aap,subj,avipath); end
            aap = aas_report_add(aap,subj,'</td></tr></table>');
        end
        
    case 'doit'        
        %% 0) Check that the templates and images we need exist!
        subj = varargin{1};
        
        % Process streams
        domain = aap.tasklist.currenttask.domain;
        [diagstream, mainstream, wbstream] = process_streams(aap);
        
        % Check precalculated XFM
        if aas_stream_has_contents(aap,'subject',subj,'epitotarget_xfm')

            load(aas_getfiles_bystream(aap,'subject',subj,'epitotarget_xfm'),'xfm');
            MMs = spm_matrix(xfm);
        end
            
        % Check local structural
        Simg = aas_getfiles_bystream(aap,subj,'structural');
        if size(Simg,1) > 1
            aas_log(aap, false, 'WARNING: Found more than 1 structural images, using first.');
            Simg = [deblank(Simg(1,:)),',1'];
        end
        % - BET structural
        Sbrain = spm_file(Simg,'suffix','_brain');
        aas_runfslcommand(aap,['bet ' Simg ' ' Sbrain ' -R']);

        % Check WM segmentation
        WMseg = '';
        if aas_stream_has_contents(aap,'subject',subj,'native_white')
            WMseg = aas_getfiles_bystream(aap,'subject',subj,'native_white');
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
        MM0 = spm_get_space(mEPIimg);

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
        load(aas_getfiles_bystream(aap,subj,'t1totemplate_xfm'),'xfm');

        aas_log(aap,false,sprintf(['\tto template realignment parameters:\n' ...
            '\tx: %0.4f   y: %0.4f   z: %0.4f   p: %0.4f   r: %0.4f   j: %0.4f'], ...
            xfm(1), xfm(2), xfm(3), xfm(4), xfm(5), xfm(6)))

        if ~exist('MMs','var')
            %% 1) Mean Functional to T1 template (reorient)
            % Set the new space for the mean functional
            for img = {mEPIimg WBimg}
                if isempty(img{1}), continue; end
                spm_get_space(img{1}, spm_matrix(xfm)\spm_get_space(img{1}));
            end
            
            %% 2) Mean Functional to Structural (coregister)
            
            if ~isempty(WBimg) % two-step coreg
                aas_log(aap,true,'NYI - epireg with whole-brain');
            else
                % Coregister mean EPI to structural
                cmd = ['epi_reg --epi=' mEPIimg ' --t1=' Simg ' --t1brain=' Sbrain ' --out=' spm_file(mEPIimg,'suffix','_reg')];
                if ~isempty(WMseg), cmd = [cmd ' --wmseg=' WMseg]; end
                aas_runfslcommand(aap,cmd);
                Mtrans = fslmat2transmat(spm_file(mEPIimg,'suffix','_reg','ext','mat'),mEPIimg,Sbrain);
                x = spm_imatrix(Mtrans);
                spm_get_space(mEPIimg,Mtrans\spm_get_space(mEPIimg));
            end
            
            aas_log(aap,false,sprintf(['\tmean EPI to structural realignment parameters:\n' ...
                '\tx: %0.4f   y: %0.4f   z: %0.4f   p: %0.4f   r: %0.4f   j: %0.4f'], ...
                x(1), x(2), x(3), x(4), x(5), x(6)))
            
            %% 3) Now apply this transformation to all the EPI images
            % The mean EPI will already be in the space required for the
            % individual EPIs. Hence, we can...
            
            MMe = spm_get_space(mEPIimg);
        end
        
        % Locate all the EPIs we want to coregister
        for m = 1:numel(mainstream)
            if ~aas_stream_has_contents(aap,domain,cell2mat(varargin),mainstream{m}), continue; end
            EPIimg{m} = cellstr(aas_getfiles_bystream(aap,domain,cell2mat(varargin),mainstream{m}));
            excl = [];
            for e = 1:numel(EPIimg{m})
                isOK = true; try spm_vol(EPIimg{m}{e}); catch, isOK = false; end
                if ~isOK
                    aas_log(aap,false,sprintf('WARNING: file %s is not a NIfTI --> skipping',EPIimg{m}{e}));
                    excl(end+1) = e;
                    continue; 
                end
                % Apply the space of the coregistered mean EPI to the
                % remaining EPIs (safest solution!)
                if exist('MMs','var')
                    MMe0 = spm_get_space(EPIimg{m}{e});
                    MMe = MMs\MMe0;
                end
                spm_get_space(EPIimg{m}{e}, MMe);
            end
            EPIimg{m}(excl) = [];
        end
        
        %% Describe the outputs and Diagnostics        
        if ~isfield(aap.tasklist.currenttask.settings,'diagnostic') || ~isstruct(aap.tasklist.currenttask.settings.diagnostic) && aap.tasklist.currenttask.settings.diagnostic
            aas_checkreg(aap,domain,cell2mat(varargin),mEPIimg,'structural');
            for m = 1:numel(mainstream)
                aas_checkreg(aap,domain,cell2mat(varargin),mainstream{m},'structural');
            end
        else % struct aap.tasklist.currenttask.settings.diagnostic
            if isfield(aap.tasklist.currenttask.settings.diagnostic,'streamind')
                aas_checkreg(aap,domain,cell2mat(varargin),mainstream{aap.tasklist.currenttask.settings.diagnostic.streamind},'structural');
            end
        end
        
        xfm = spm_imatrix(inv(MMe/MM0));
        fn = fullfile(aas_getpath_bydomain(aap,domain,cell2mat(varargin)),'epitotarget_xfm.mat');
        save(fn,'xfm');
        aap = aas_desc_outputs(aap,domain,cell2mat(varargin),'epitotarget_xfm',fn);
        
        for m = 1:numel(mainstream)
            if ~aas_stream_has_contents(aap,domain,cell2mat(varargin),mainstream{m}), continue; end
            aap = aas_desc_outputs(aap,domain,cell2mat(varargin),mainstream{m},EPIimg{m});
        end
        if any(strcmp(aas_getstreams(aap,'output'),diagstream))
            aap = aas_desc_outputs(aap,domain,cell2mat(varargin),diagstream,mEPIimg);
        end
        
    case 'checkrequirements'
        [in, attr] = aas_getstreams(aap,'input'); 
        indModRef = find(cellfun(@(a) isstruct(a) && isfield(a,'diagnostic') && a.diagnostic, attr),1,'last');
        in(1:indModRef-1) = []; % not for cross-modality reference(s)
        out = aas_getstreams(aap,'output'); if ~iscell(out), out = {out}; end 
        out = setdiff(out,{'epitotarget_xfm'},'stable');
        for s = 1:numel(in)
            instream = strsplit(in{s},'.'); instream = instream{end};
            if s <= numel(out)
                if ~strcmp(out{s},instream)
                    if s == 1 % modality reference missing (optional)
                        out = [{''} out];
                        continue
                    end
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
outpstreams = aas_getstreams(aap,'output'); outpstreams = setdiff(outpstreams,'epitotarget_xfm');

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
[~,ia] = intersect(inpstreams,outpstreams);
outind = false(size(mainind)); outind(ia) = true;
mainind = mainind & outind;
if any(mainind)
    mainstream = inpstreams(mainind);
else
    mainstream = {};
end
end

function transmat = fslmat2transmat(fslmat, src, trg)
% This is a slightly modified version of flirtmat2worldmat
% by Ged Ridgway available at 
% https://www.nitrc.org/snippet/detail.php?type=snippet&id=5

src = spm_vol(src);
trg = spm_vol(trg);
fslmat = dlmread(fslmat);

% src.mat = flirtmat \ trg.mat
% srcvox = src.mat \ inv(flirtmat) * trg.mat * trgvox
% BUT, flirt doesn't use src.mat, only absolute values of the 
% scaling elements from it,
% AND, if images are not radiological, the x-axis is flipped, see:
%  https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind0810&L=FSL&P=185638
%  https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind0903&L=FSL&P=R93775
fslvoxmat = nifti2scl(src) \ inv(fslmat) * nifti2scl(trg);

% AND, Flirt's voxels are zero-based, while SPM's are one-based...
addone = eye(4); addone(:, 4) = 1;
spmvoxmat = addone * fslvoxmat / addone;

transmat = src.mat * spmvoxmat / trg.mat;
end

%%
function scl = nifti2scl(V)
% not sure if this is always correct with rotations in mat, but seems okay!
scl = diag([sqrt(sum(V.mat(1:3,1:3).^2)) 1]);
if det(V.mat) > 0
    % neurological, x-axis is flipped, such that [3 2 1 0] and [0 1 2 3]
    % have the same *scaled* coordinates:
    xflip = diag([-1 1 1 1]);
    xflip(1, 4) = V.dim(1)-1; % reflect about centre
    scl = scl * xflip;
end
end
