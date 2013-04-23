function cg_vbm8_defaults
% Sets the defaults for VBM
% FORMAT cg_vbm_defaults
%_______________________________________________________________________
%
% This file is intended to be customised for the site.
%
% Care must be taken when modifying this file
%_______________________________________________________________________
% $Id: cg_vbm8_defaults.m 432 2011-10-03 21:48:53Z gaser $

global vbm8

% Estimation options
%=======================================================================
vbm8.opts.tpm       = {fullfile(spm('dir'),'toolbox','Seg','TPM.nii')}; % TPM.nii
vbm8.opts.ngaus     = [2 2 2 3 4 2];  % Gaussians per class
vbm8.opts.affreg    = 'mni';    % Affine regularisation
vbm8.opts.warpreg   = 4;      % Warping regularisation
vbm8.opts.biasreg   = 0.0001; % Bias regularisation
vbm8.opts.biasfwhm  = 60;   % Bias FWHM
vbm8.opts.samp      = 3;      % Sampling distance

% Writing options
%=======================================================================

% segmentations:
%   native    0/1   (none/yes)
%   warped    0/1   (none/yes)
%   modulated 0/1/2 (none/affine+nonlinear/nonlinear only)
%   dartel    0/1/2 (none/rigid/affine)

vbm8.output.bias.native  = 0;
vbm8.output.bias.warped  = 1;
vbm8.output.bias.affine  = 0;

vbm8.output.label.native = 0;
vbm8.output.label.warped = 0;
vbm8.output.label.dartel = 0;

% order is [native normalised modulated dartel]
vbm8.output.GM.native = 0;  % GM
vbm8.output.GM.warped = 0;  % GM
vbm8.output.GM.mod    = 2;  % GM
vbm8.output.GM.dartel = 0;  % GM

vbm8.output.WM.native = 0;  % WM
vbm8.output.WM.warped = 0;  % WM
vbm8.output.WM.mod    = 2;  % WM
vbm8.output.WM.dartel = 0;  % WM

vbm8.output.CSF.native = 0; % CSF
vbm8.output.CSF.warped = 0; % CSF
vbm8.output.CSF.mod    = 0; % CSF
vbm8.output.CSF.dartel = 0; % CSF

% jacobian determinant 0/1 (none/yes)
vbm8.output.jacobian.warped = 0;

% order is [forward inverse]
vbm8.output.warps = [0 0];

% Extended writing options
%=======================================================================
vbm8.extopts.dartelwarp  = 1; % dartel normalization: 0 - spm default; 1 - yes
vbm8.extopts.darteltpm   = {fullfile(spm('dir'),'toolbox','vbm8','Template_1_IXI550_MNI152.nii')}; % Indicate first Dartel template
vbm8.extopts.print       = 1; % Display and print results

% bias correction options
%=======================================================================
vbm8.bias.nits_bias    = 8;
vbm8.bias.biasfwhm     = 60;
vbm8.bias.biasreg      = 1e-6;
vbm8.bias.lmreg        = 1e-6;

% realign options
%=======================================================================
vbm8.realign.halfway   = 1; % use halfway registration: 0 - no; 1 - yes
vbm8.realign.weight    = 1; % weight registration with inverse std: 0 - no; 1 - yes
vbm8.realign.ignore_mat= 0; % ignore exisiting positional information: 0 - no; 1 - yes

% apply deformations options
%=======================================================================
vbm8.defs.interp    = 5;  % 5th degree B-spline

% expert options
%=======================================================================
vbm8.extopts.cleanup   = 1;    % Cleanup: 0 - no; 1 - light; 2 -thorough
vbm8.extopts.finalmask = 1;    % Final masking: 0 - no; 1 - yes
vbm8.extopts.gcut      = 1;    % Skull-stripping with graph-cut: 0 - no; 1 - yes
vbm8.extopts.kmeans    = 1;    % segmentation initialization: 0 - new segment; 1 - Kmeans
vbm8.extopts.mrf       = 0.15; % MRF weighting
vbm8.extopts.sanlm     = 2;    % use SANLM filter: 0 - no SANLM; 1 - SANLM with single-threading; 2 - SANLM with multi-threading
vbm8.extopts.bias_fwhm = 60;   % FWHM of Kmeans internal bias correction

% experimental (not yet working)
%=======================================================================
vbm8.extopts.mask      = {fullfile(spm('dir'),'toolbox','vbm8','submask.nii')}; % mask for subcortical areas + ventricles
vbm8.output.surf.dartel= 0; % WM-surface 
