function tfce_estimate = tbx_cfg_tfce_estimate
% SPM Configuration file for TFCE estimate
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: tbx_cfg_tfce_estimate.m 82 2015-11-04 13:07:20Z gaser $

rev = '$Rev: 82 $';

addpath(fileparts(which(mfilename)));

% ---------------------------------------------------------------------
% spmmat Select SPM.mat
% ---------------------------------------------------------------------
spmmat         = cfg_files;
spmmat.tag     = 'spmmat';
spmmat.name    = 'Select SPM.mat';
spmmat.help    = {'Select the SPM.mat file that contains the design specification.'};
spmmat.filter  = 'mat';
spmmat.ufilter = '^SPM\.mat$';
spmmat.num     = [1 1];

% ---------------------------------------------------------------------
% mask Select mask image to restrict analysis
% ---------------------------------------------------------------------
mask         = cfg_files;
mask.tag     = 'mask';
mask.name    = 'Select additional mask image';
mask.help    = {'Select an additional mask image to restrict analysis. As default the mask image in the analysis folder is used. Here you can select a mask image to additionally restrict the analysis to regions of interest (i.e. small volume correction).'};
if strcmp(spm('ver'),'SPM12')
  mask.filter  = {'image','mesh'};
else
  mask.filter  = {'image'};
end
mask.val     = {''};
mask.ufilter = '.*';
mask.num     = [0 1];

% ---------------------------------------------------------------------
% titlestr Results Title
% ---------------------------------------------------------------------
titlestr         = cfg_entry;
titlestr.tag     = 'titlestr';
titlestr.name    = 'Results Title';
titlestr.help    = {'Heading on results page - determined automatically if left empty'};
titlestr.val     = {''};
titlestr.strtype = 's';
titlestr.num     = [0 Inf];

% ---------------------------------------------------------------------
% contrasts Contrast
% ---------------------------------------------------------------------
contrasts         = cfg_entry;
contrasts.tag     = 'contrasts';
contrasts.name    = 'Contrast index';
contrasts.help    = {'Index(es) of contrast according to the contrast manager.'
                     ''
                     'Each contrast in SPM is indicated by a sequential number that is displayed in the first column of the contrast manager.'
                     ''
                     'You can enter one or more contrasts. If only one number is entered, and this number is "Inf", you can select one or more contrasts interactively using the contrast manager.'
                     ''
                     'Do not define here the contrast itself. This should be done in the contrast manager, that is automatically called if "Inf" is kept as entry.'
}';
contrasts.strtype = 'e';
contrasts.val     = {Inf};
contrasts.num     = [1 Inf];

% ---------------------------------------------------------------------
% number of permutations
% ---------------------------------------------------------------------
n_perm         = cfg_entry;
n_perm.tag     = 'n_perm';
n_perm.name    = 'Number of permutations';
n_perm.help    = {'With 1000 permutations the smallest possible p-value is 0.001 (n=1/p). A useful strategy is to start with 1000 permutations and continue to 5000-10000 only if p is small enough to be interesting and/or for the final analysis.'
                     ''
                     'If number of maximal possible permutations is smaller, then this number is used resulting in an exact permutation test.'
}';
n_perm.strtype = 'e';
n_perm.val     = {5000};
n_perm.num     = [1 Inf];

% ---------------------------------------------------------------------
% two-dimensional processing
% ---------------------------------------------------------------------
tbss    = cfg_menu;
tbss.tag = 'tbss';
tbss.name = 'TBSS data';
tbss.labels = {'yes','no'};
tbss.values = {1 0};
tbss.val  = {0};
tbss.help = {[...
'Use 2D optimization (e.g. for TBSS data) with internal TFCE parameters H=2, E=1.']};

% ---------------------------------------------------------------------
% variance smoothing
% ---------------------------------------------------------------------
vFWHM         = cfg_entry;
vFWHM.tag     = 'vFWHM';
vFWHM.name    = 'Variance smoothing (for low DFs)';
vFWHM.help    = {'For low degrees of freedom (n<20) the variance is estimated poorly and smoothing the variance can help in obtaining a more reliable estimate. By smoothing the variance the noise from the variance will be smoothed, but not the signal which results in a so-called pseudo t-statistics.'
}';
vFWHM.strtype = 'e';
vFWHM.val     = {0};
vFWHM.num     = [0 Inf];

% ---------------------------------------------------------------------
% conspec Contrast query
% ---------------------------------------------------------------------
conspec         = cfg_branch;
conspec.tag     = 'conspec';
conspec.name    = 'Contrast query';
conspec.val     = {titlestr contrasts n_perm vFWHM};
conspec.help    = {''};

% ---------------------------------------------------------------------
% multithreading
% ---------------------------------------------------------------------
openmp    = cfg_menu;
openmp.tag = 'openmp';
openmp.name = 'Use multi-threading to speed up calculations';
openmp.labels = {'yes','no'};
openmp.values = {1 0};
openmp.val  = {1};
openmp.help = {[...
'OpenMP can be used to distribute calculations to multiple processors. ',...
'This will minimize calculation time by a large amount, but makes sometimes trouble on Windows machines. ',...
'In case of trouble deselect the multi-threading option.']};

% ---------------------------------------------------------------------
% results Results Report
% ---------------------------------------------------------------------
tfce_estimate          = cfg_exbranch;
tfce_estimate.tag      = 'tfce_estimate';
tfce_estimate.name     = 'Estimate TFCE';
tfce_estimate.val      = {spmmat mask conspec tbss openmp};
tfce_estimate.help     = {''};
tfce_estimate.prog     = @cg_tfce_estimate;
