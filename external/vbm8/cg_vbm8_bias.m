function bias = cg_vbm8_bias
% Configuration file for bias correction between an image pair
%
% Christian Gaser
% $Id: cg_vbm8_bias.m 404 2011-04-11 10:03:40Z gaser $

mov = cfg_files;
mov.name = 'Longitudinal images for one subject';
mov.tag  = 'mov';
mov.filter = 'image';
mov.num  = [1 Inf];
mov.help   = {[...
'These are the images of the same subject. The first image is used as reference and all subsequent images are bias corrected with regard to the first image']};
%------------------------------------------------------------------------

subj = cfg_branch;
subj.name = 'Subject';
subj.tag = 'subj';
subj.val = {mov};
subj.help = {[...
'Two images of the same subject, which are to be registered together.  Prior to bias correction, the images should be rigidly registered with each other.']};

%------------------------------------------------------------------------

esubjs         = cfg_repeat;
esubjs.tag     = 'esubjs';
esubjs.name    = 'Data';
esubjs.values  = {subj };
esubjs.num     = [1 Inf];
esubjs.help = {[...
'Specify pairs of images to match together.']};

%------------------------------------------------------------------------

nits = cfg_entry;
nits.name = 'Iterations for bias correction';
nits.tag  = 'nits';
nits.strtype = 'n';
nits.num  = [1 1];
nits.def     = @(val)cg_vbm8_get_defaults('bias.nits_bias', val{:});
nits.help = {'Number of iterations for the bias correction.'};

%------------------------------------------------------------------------

biasfwhm = cfg_menu;
biasfwhm.name = 'Bias FWHM';
biasfwhm.tag  = 'fwhm';
biasfwhm.labels = {...
'30mm cutoff','40mm cutoff','50mm cutoff','60mm cutoff','70mm cutoff',...
'80mm cutoff','90mm cutoff','100mm cutoff','110mm cutoff','120mm cutoff',...
'130mm cutoff','140mm cutoff','150mm cutoff','No correction'};
biasfwhm.values = {30,40,50,60,70,80,90,100,110,120,130,140,150,Inf};
biasfwhm.def     = @(val)cg_vbm8_get_defaults('bias.biasfwhm', val{:});
biasfwhm.help = {[...
'FWHM of Gaussian smoothness of bias. If your intensity nonuniformity is very smooth, then choose a large FWHM. This will prevent the algorithm from trying to model out intensity variation due to different tissue types. The model for intensity nonuniformity is one of i.i.d. Gaussian noise that has been smoothed by some amount, before taking the exponential. Note also that smoother bias fields need fewer parameters to describe them. This means that the algorithm is faster for smoother intensity nonuniformities.']};

%------------------------------------------------------------------------

biasreg = cfg_menu;
biasreg.name = 'Bias regularisation';
biasreg.tag = 'reg';
biasreg.labels = {...
'no regularisation','extremely light regularisation',...
'very light regularisation','light regularisation',...
'medium regularisation','heavy regularisation',...
'very heavy regularisation','extremely heavy regularisation'};
biasreg.values = {0, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3};
biasreg.def     = @(val)cg_vbm8_get_defaults('bias.biasreg', val{:});
biasreg.help = {[...
'We know a priori that intensity variations due to MR physics tend to be spatially smooth, whereas those due to different tissue types tend to contain more high frequency information. A more accurate estimate of a bias field can be obtained by including prior knowledge about the distribution of the fields likely to be encountered by the correction algorithm. For example, if it is known that there is little or no intensity non-uniformity, then it would be wise to penalise large values for the intensity nonuniformity parameters. This regularisation can be placed within a Bayesian context, whereby the penalty incurred is the negative logarithm of a prior probability for any particular pattern of nonuniformity.']};

%------------------------------------------------------------------------

lmreg = cfg_menu;
lmreg.name = 'Levenberg-Marquardt regularisation';
lmreg.tag = 'lmreg';
lmreg.labels = {...
'no regularisation','extremely light regularisation',...
'very light regularisation','light regularisation',...
'medium regularisation','heavy regularisation',...
'very heavy regularisation','extremely heavy regularisation'};
lmreg.values = {0, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3};
lmreg.def     = @(val)cg_vbm8_get_defaults('bias.lmreg', val{:});
lmreg.help = {[...
'Levenberg-Marquardt regularisation keeps the bias correction part stable. Higher values means more stability, but slower convergence.']};

%------------------------------------------------------------------------

bias_opts = cfg_branch;
bias_opts.name = 'Bias Correction Options';
bias_opts.tag = 'bias_opts';
bias_opts.val = {nits,biasfwhm,biasreg,lmreg};
bias_opts.help = {[...
'MR images are usually corrupted by a smooth, spatially varying artifact that modulates the intensity of the image (bias). These artifacts, although not usually a problem for visual inspection, can impede automated processing of the images.'],...
[...
'Before registering the images, an approximate bias correction is estimated for the moved image. This is based on minimising the difference between the images an a symmetric way. Prior to registering the images, they should be rigidly aligned together.  The bias correction is estimated once for these aligned images.']};

%------------------------------------------------------------------------

bias = cfg_exbranch;
bias.name = 'Intra-subject bias correction';
bias.tag  = 'bias';
bias.val  = {esubjs bias_opts};
bias.prog = @cg_vbm8_bias_run;
bias.vout = @vout_bias;
bias.help = {
'This option provides an intra-subject bias correction. The first image is used as reference to correct bias of all subsequent images of the same subject.'};

%------------------------------------------------------------------------

return;
%------------------------------------------------------------------------
 
%------------------------------------------------------------------------
function dep = vout_bias(job)
for k=1:numel(job.subj)
    dep(k)            = cfg_dep;
    dep(k).sname      = sprintf('Bias corrected images (Subj %d)',k);
    dep(k).src_output = substruct('()',{k},'.','files');
    dep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end;
%------------------------------------------------------------------------
