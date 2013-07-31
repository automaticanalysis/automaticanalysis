function varargout = cg_vbm8_run(job,arg)
% Segment a bunch of images
% FORMAT cg_vbm8_run(job)
% job.channel(n).vols{m}
% job.channel(n).biasreg
% job.channel(n).biasfwhm
% job.channel(n).write
% job.tissue(k).tpm
% job.tissue(k).ngaus
% job.tissue(k).native
% job.tissue(k).warped
% job.warp.affreg
% job.warp.reg
% job.warp.samp
% job.warp.write
% job.warp.dartelwarp
% job.warp.print
%
% See the user interface for a description of the fields.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% based on John Ashburners version of
% spm_preproc8_run.m 2281 2008-10-01 12:52:50Z john $
%
% Christian Gaser
% $Id: cg_vbm8_run.m 411 2011-04-14 13:39:03Z gaser $

rev = '$Rev: 411 $';

% check whether estimation & write should be done
estwrite = isfield(job,'opts');

% set some defaults if segmentations are not estimated
if ~estwrite
    job.opts = struct('biasreg',0.001,'biasfwhm',60,'affreg','mni',...
                      'warpreg',4,'samp',3,'ngaus',[2 2 2 3 4 2]);
end

channel = struct('vols',{job.data});
                 
warp = struct('affreg', job.opts.affreg,...
              'samp', job.opts.samp,...
              'reg', job.opts.warpreg,...
              'write', job.output.warps,...
              'sanlm', job.extopts.sanlm,...
              'mrf', job.extopts.mrf,...
              'print', job.extopts.print,...
              'cleanup', job.extopts.cleanup,...
              'dartelwarp', isfield(job.extopts.dartelwarp,'normhigh'));

if isfield(job.extopts.dartelwarp,'normhigh')
    warp.darteltpm = job.extopts.dartelwarp.normhigh.darteltpm{1};
end

do_dartel = warp.dartelwarp;

% prepare tissue priors and number of gaussians for all 6 classes
if estwrite
    [pth,nam,ext,num] = spm_fileparts(job.opts.tpm{1});
    for i=1:6
        tissue(i).ngaus = job.opts.ngaus(i);
        tissue(i).tpm = [fullfile(pth,[nam ext]) ',' num2str(i)];
    end
end

% write tissue class 1-3              
tissue(1).warped = [job.output.GM.warped  (job.output.GM.modulated==1)  (job.output.GM.modulated==2) ];
tissue(1).native = [job.output.GM.native  (job.output.GM.dartel==1)     (job.output.GM.dartel==2)    ];
tissue(2).warped = [job.output.WM.warped  (job.output.WM.modulated==1)  (job.output.WM.modulated==2) ];
tissue(2).native = [job.output.WM.native  (job.output.WM.dartel==1)     (job.output.WM.dartel==2)    ];
tissue(3).warped = [job.output.CSF.warped (job.output.CSF.modulated==1) (job.output.CSF.modulated==2)];
tissue(3).native = [job.output.CSF.native (job.output.CSF.dartel==1)    (job.output.CSF.dartel==2)   ];

% never write class 4-6
for i=4:6
    tissue(i).warped = [0 0 0];
    tissue(i).native = [0 0 0];
end

job.bias     = [job.output.bias.native  job.output.bias.warped job.output.bias.affine];
job.label    = [job.output.label.native job.output.label.warped (job.output.label.dartel==1) (job.output.label.dartel==2)];
job.jacobian = job.output.jacobian.warped;

job.biasreg  = job.opts.biasreg;
job.biasfwhm = job.opts.biasfwhm;
job.channel  = channel;
job.warp     = warp;
job.warps    = job.output.warps;
job.tissue   = tissue;

if nargin == 1, arg = 'run'; end

switch lower(arg)
    case 'run'
        varargout{1} = run_job(job,estwrite);
    case 'check'
        varargout{1} = check_job(job);
    case 'vfiles'
        varargout{1} = vfiles_job(job);
    case 'vout'
        varargout{1} = vout_job(job);
    otherwise
        error('Unknown argument ("%s").', arg);
end
return
%_______________________________________________________________________

%_______________________________________________________________________
function vout = run_job(job, estwrite)

vout   = vout_job(job);

% load tpm priors only for estimate and write
if estwrite
    tpm    = strvcat(cat(1,job.tissue(:).tpm));
    tpm    = spm_load_priors8(tpm);
end

nit = 1;

for iter=1:nit,
    if nit>1,
        % Sufficient statistics for possible generation of group-specific
        % template data.
        SS = zeros([size(tpm.dat{1}),numel(tpm.dat)],'single');
    end
    for subj=1:numel(job.channel(1).vols),
        if estwrite % estimate and write segmentations
            images = '';
            for n=1:numel(job.channel),
                images = strvcat(images,job.channel(n).vols{subj});
            end
            obj.image    = spm_vol(images);
            spm_check_orientations(obj.image);

            obj.fudge    = 5;
            obj.biasreg  = cat(1,job.biasreg);
            obj.biasfwhm = cat(1,job.biasfwhm);
            obj.tpm      = tpm;
            obj.lkp      = [];
            if all(isfinite(cat(1,job.tissue.ngaus))),
                for k=1:numel(job.tissue),
                    obj.lkp = [obj.lkp ones(1,job.tissue(k).ngaus)*k];
                end;
            end
            obj.reg      = job.warp.reg;
            obj.samp     = job.warp.samp;

            if iter==1,
                % Initial affine registration.
                Affine  = eye(4);
                if ~isempty(job.warp.affreg),
                    VG = spm_vol(fullfile(spm('Dir'),'templates','T1.nii'));
                    VF = spm_vol(obj.image(1));
                    
                    % smooth source with 8mm
                    VF1 = spm_smoothto8bit(VF,8);

                    % Rescale images so that globals are better conditioned
                    VF1.pinfo(1:2,:) = VF1.pinfo(1:2,:)/spm_global(VF1);
                    VG.pinfo(1:2,:)  = VG.pinfo(1:2,:)/spm_global(VG);

                    fprintf('Initial Coarse Affine Registration..\n');
                    aflags    = struct('sep',8, 'regtype',job.warp.affreg,...
                        'WG',[],'WF',[],'globnorm',0);
                    aflags.sep = max(aflags.sep,max(sqrt(sum(VG(1).mat(1:3,1:3).^2))));
                    aflags.sep = max(aflags.sep,max(sqrt(sum(VF(1).mat(1:3,1:3).^2))));

                    M = eye(4);
                    try
                        spm_plot_convergence('Init','Coarse Affine Registration','Mean squared difference','Iteration');
                    catch
                        spm_chi2_plot('Init','Coarse Affine Registration','Mean squared difference','Iteration');
                    end
                    [Affine, scale]  = spm_affreg(VG, VF1, aflags, M);

                    aflags.WG  = spm_vol(fullfile(spm('Dir'),'apriori','brainmask.nii'));
                    aflags.sep = aflags.sep/2;
                    try
                        spm_plot_convergence('Init','Fine Affine Registration','Mean squared difference','Iteration');
                    catch
                        spm_chi2_plot('Init','Fine Affine Registration','Mean squared difference','Iteration');
                    end
                    Affine  = spm_affreg(VG, VF1, aflags, Affine, scale);

                    fprintf('Fine Affine Registration..\n');
                    Affine  = spm_maff8(obj.image(1),job.warp.samp,obj.fudge,  tpm,Affine,job.warp.affreg);
                end;
                obj.Affine = Affine;
            else
                % Load results from previous iteration for use with next round of
                % iterations, with the new group-specific tissue probability map.
                [pth,nam] = fileparts(job.channel(1).vols{subj});
                res       = load(fullfile(pth,[nam '_seg8.mat']));
                obj.Affine = res.Affine;
                obj.Twarp  = res.Twarp;
                obj.Tbias  = res.Tbias;
                if ~isempty(obj.lkp),
                    obj.mg     = res.mg;
                    obj.mn     = res.mn;
                    obj.vr     = res.vr;
                end
            end

            res = spm_preproc8(obj);

            try
                [pth,nam] = fileparts(job.channel(1).vols{subj});
                savefields(fullfile(pth,[nam '_seg8.mat']),res);
            catch
            end

        else % only write segmentations
            [pth,nam] = fileparts(job.channel(1).vols{subj});
            seg8_name = fullfile(pth,[nam '_seg8.mat']);
            if exist(seg8_name)
                res = load(seg8_name);

                % load original used tpm, which is save in seg8.mat file
                try
                    tpm    = spm_load_priors8(res.tpm);
                catch
                    % or use default TPM
                    fprintf('Original TPM image %s was not found. use default TPM image instead.\n',res.tpm(1).fname);
                    for i=1:6
                      job.tissue(i).tpm = fullfile(spm('dir'),'toolbox','Seg',['TPM.nii,' num2str(i)]);
                    end
                    tpm    = strvcat(cat(1,job.tissue(:).tpm));
                    tpm    = spm_load_priors8(tpm);
                end
                
                % use path of mat-file in case that image was moved
				        [image_pth,image_nam,image_ext] = spm_fileparts(job.channel(1).vols{subj});
				        res.image(1).fname = fullfile(image_pth, [image_nam, image_ext]);
            else
                error(['Can''t load file ' seg8_name]);  
                return
            end
        end
        if iter==nit,
            % Final iteration, so write out the required data.
            tc = [cat(1,job.tissue(:).native) cat(1,job.tissue(:).warped)];
            bf = job.bias;
            df = job.warp.write;
            lb = job.label;
            jc = job.jacobian;
            cg_vbm8_write(res, tc, bf, df, lb, jc, job.warp, tpm, job)
        else
            % Not the final iteration, so compute sufficient statistics for
            % re-estimating the template data.
            N    = numel(job.channel);
            K    = numel(job.tissue);
            cls  = cg_vbm8_write(res,zeros(K,4),zeros(N,2),[0 0],[0 0 0 0], 0, job.warp, tpm, job);
            for k=1:K,
                SS(:,:,:,k) = SS(:,:,:,k) + cls{k};
            end
        end

    end
    if iter<nit && nit>1,
         % Treat the tissue probability maps as Dirichlet priors, and compute the 
         % MAP estimate of group tissue probability map using the sufficient
         % statistics.
         alpha = 1.0;
         for k=1:K,
             SS(:,:,:,k) = SS(:,:,:,k) + spm_bsplinc(tpm.V(k),[0 0 0  0 0 0])*alpha + eps;
         end

         s = sum(SS,4);
         for k=1:K,
             tmp        = SS(:,:,:,k)./s;
             tpm.bg1(k) = mean(mean(tmp(:,:,1)));
             tpm.bg2(k) = mean(mean(tmp(:,:,end)));
             tpm.dat{k} = spm_bsplinc(log(tmp+tpm.tiny),[ones(1,3)*(tpm.deg-1)  0 0 0]);
         end
    end
end
return
%_______________________________________________________________________

%_______________________________________________________________________
function msg = check_job(job)
msg = {};
if numel(job.channel) >1,
    k = numel(job.channel(1).vols);
    for i=2:numel(job.channel),
        if numel(job.channel(i).vols)~=k,
            msg = {['Incompatible number of images in channel ' num2str(i)]};
            break
        end
    end
elseif numel(job.channel)==0,
    msg = {'No data'};
end
return
%_______________________________________________________________________

%_______________________________________________________________________
function savefields(fnam,p)
if length(p)>1, error('Can''t save fields.'); end;
fn = fieldnames(p);
if numel(fn)==0, return; end;
for i=1:length(fn),
    eval([fn{i} '= p.' fn{i} ';']);
end;
try 
    mat_ver = spm_check_version('matlab','7');
catch
    mat_ver = spm_matlab_version_chk('7');
end 
if mat_ver >= 0
    save(fnam,'-V6',fn{:});
else
    save(fnam,fn{:});
end;

return;
%_______________________________________________________________________

%_______________________________________________________________________
function vout = vout_job(job)

do_dartel = isfield(job.extopts.dartelwarp,'normhigh');

n     = numel(job.channel(1).vols);
parts = cell(n,4);

biascorr  = {};
wbiascorr = {};
label  = {};
wlabel = {};
rlabel = {};
alabel = {};
jacobian = {};

for j=1:n,
    [parts{j,:}] = spm_fileparts(job.channel(1).vols{j});
end

if job.bias(1),
    biascorr = cell(n,1);
    for j=1:n
        biascorr{j} = fullfile(parts{j,1},['m',parts{j,2},'.nii']);
    end
end

if job.bias(2),
    wbiascorr = cell(n,1);
    for j=1:n
        if do_dartel
            wbiascorr{j} = fullfile(parts{j,1},['wmr',parts{j,2},'.nii']);
        else
            wbiascorr{j} = fullfile(parts{j,1},['wm',parts{j,2},'.nii']);
        end
    end
end

if job.label(1),
    label = cell(n,1);
    for j=1:n
        label{j} = fullfile(parts{j,1},['p0',parts{j,2},'.nii']);
    end
end

if job.label(2),
    wlabel = cell(n,1);
    for j=1:n
        if do_dartel
            wlabel{j} = fullfile(parts{j,1},['wrp0',parts{j,2},'.nii']);
        else
            wlabel{j} = fullfile(parts{j,1},['wp0',parts{j,2},'.nii']);
        end
    end
end

if job.label(3),
    rlabel = cell(n,1);
    for j=1:n
        rlabel{j} = fullfile(parts{j,1},['rp0',parts{j,2},'.nii']);
    end
end

if job.label(4),
    alabel = cell(n,1);
    for j=1:n
        alabel{j} = fullfile(parts{j,1},['rp0',parts{j,2},'_affine.nii']);
    end
end

for j=1:n,
    [parts{j,:}] = spm_fileparts(job.channel(1).vols{j});
end
param = cell(n,1);
for j=1:n
    param{j} = fullfile(parts{j,1},[parts{j,2},'_seg8.mat']);
end

tiss = struct('c',{},'rc',{},'rca',{},'wc',{},'mwc',{},'m0wc',{});
for i=1:numel(job.tissue),
    if job.tissue(i).native(1),
        tiss(i).c = cell(n,1);
        for j=1:n
            tiss(i).c{j} = fullfile(parts{j,1},['p',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).native(2),
        tiss(i).rc = cell(n,1);
        for j=1:n
            tiss(i).rc{j} = fullfile(parts{j,1},['rp',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).native(3),
        tiss(i).rca = cell(n,1);
        for j=1:n
            tiss(i).rca{j} = fullfile(parts{j,1},['rp',num2str(i),parts{j,2},'_affine.nii']);
        end
    end
    if job.tissue(i).warped(1),
        tiss(i).wc = cell(n,1);
        for j=1:n
            if do_dartel
                tiss(i).wc{j} = fullfile(parts{j,1},['wrp',num2str(i),parts{j,2},'.nii']);
            else
                tiss(i).wc{j} = fullfile(parts{j,1},['wp',num2str(i),parts{j,2},'.nii']);
            end
        end
    end
    if job.tissue(i).warped(2),
        tiss(i).mwc = cell(n,1);
        for j=1:n
            if do_dartel
                tiss(i).mwc{j} = fullfile(parts{j,1},['mwrp',num2str(i),parts{j,2},'.nii']);
            else
                tiss(i).mwc{j} = fullfile(parts{j,1},['mwp',num2str(i),parts{j,2},'.nii']);
            end
        end
    end
    if job.tissue(i).warped(3),
        tiss(i).m0wc = cell(n,1);
        for j=1:n
            if do_dartel
                tiss(i).m0wc{j} = fullfile(parts{j,1},['m0wrp',num2str(i),parts{j,2},'.nii']);
            else
                tiss(i).m0wc{j} = fullfile(parts{j,1},['m0wp',num2str(i),parts{j,2},'.nii']);
            end
        end
    end
end

if job.warp.write(1),
    fordef = cell(n,1);
    for j=1:n
        if do_dartel
            fordef{j} = fullfile(parts{j,1},['y_r',parts{j,2},'.nii']);
        else
            fordef{j} = fullfile(parts{j,1},['y_',parts{j,2},'.nii']);
        end
    end
else
    fordef = {};
end

if job.warp.write(2),
    invdef = cell(n,1);
    for j=1:n
        if do_dartel
            invdef{j} = fullfile(parts{j,1},['iy_r',parts{j,2},'.nii']);
        else
            invdef{j} = fullfile(parts{j,1},['iy_',parts{j,2},'.nii']);
        end
    end
else
    invdef = {};
end

if job.jacobian,
    jacobian = cell(n,1);
    for j=1:n
        if do_dartel
            jacobian{j} = fullfile(parts{j,1},['jac_wrp1',parts{j,2},'.nii']);
        else
            jacobian{j} = '';
        end
    end
else
    jacobian = {};
end

vout  = struct('tiss',tiss,'label',{label},'wlabel',{wlabel},'rlabel',{rlabel},'alabel',{alabel},...
               'biascorr',{biascorr},'wbiascorr',{wbiascorr},'param',{param},...
               'invdef',{invdef},'fordef',{fordef},'jacobian',{jacobian});
%_______________________________________________________________________

%_______________________________________________________________________
function vf = vfiles_job(job)
vout = vout_job(job);
vf   = vout.param;
if ~isempty(vout.invdef),     vf = {vf{:}, vout.invdef{:}}; end
if ~isempty(vout.fordef),     vf = {vf{:}, vout.fordef{:}}; end
if ~isempty(vout.jacobian),   vf = {vf{:}, vout.jacobian{:}}; end

if ~isempty(vout.biascorr),   vf = {vf{:}, vout.biascorr{:}}; end
if ~isempty(vout.wbiascorr),  vf = {vf{:}, vout.wbiascorr{:}}; end
if ~isempty(vout.label),      vf = {vf{:}, vout.label{:}}; end
if ~isempty(vout.wlabel),     vf = {vf{:}, vout.wlabel{:}}; end
if ~isempty(vout.rlabel),     vf = {vf{:}, vout.rlabel{:}}; end
if ~isempty(vout.alabel),     vf = {vf{:}, vout.alabel{:}}; end

for i=1:numel(vout.tiss)
    if ~isempty(vout.tiss(i).c),   vf = {vf{:}, vout.tiss(i).c{:}};   end
    if ~isempty(vout.tiss(i).rc),  vf = {vf{:}, vout.tiss(i).rc{:}};  end
    if ~isempty(vout.tiss(i).rca), vf = {vf{:}, vout.tiss(i).rca{:}}; end
    if ~isempty(vout.tiss(i).wc),  vf = {vf{:}, vout.tiss(i).wc{:}};  end
    if ~isempty(vout.tiss(i).mwc), vf = {vf{:}, vout.tiss(i).mwc{:}}; end
    if ~isempty(vout.tiss(i).m0wc),vf = {vf{:}, vout.tiss(i).m0wc{:}};end
end
vf = reshape(vf,numel(vf),1);
%_______________________________________________________________________

%_______________________________________________________________________

%=======================================================================
function t = transf(B1,B2,B3,T)
if ~isempty(T)
    d2 = [size(T) 1];
    t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
    t  = B1*t1*B2';
else
    t = zeros(size(B1,1),size(B2,1),size(B3,1));
end;
return;
%=======================================================================
