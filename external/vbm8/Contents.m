% Voxel Based Morphometry Toolbox
% Version  435  (VBM8)  2011-12-06
% __________________________________________________________________________
% Copyright (C) 2009 Christian Gaser christian.gaser@uni-jena.de
%
% $Id: Contents_info.txt 420 2011-07-04 14:26:39Z gaser $
% ==========================================================================
% Description
% ==========================================================================
% This toolbox is a collection of extensions to the segmentation algorithm 
% of SPM8 (Wellcome Department of Cognitive Neurology) to provide voxel-
% based morphometry (VBM). It is developed by Christian Gaser (University of 
% Jena, Department of Psychiatry) and is available to the scientific 
% community under the terms of the GNU General Public License.
%
% General files
%   INSTALL.txt                  - installation instructions
%   vbm8.man                     - notes on VBM8 toolbox
%   CHANGES.txt                  - changes in revisions
%   Contents.m                   - this file
%
% VBM8 functions
%   cg_vbm8_bias.m               - configuration file for bias correction between an image pair
%   cg_vbm8_bias_run.m           - correct bias between an image pair
%   cg_vbm8_debug.m              - print debug information for SPM8 and VBM8
%   cg_vbm8_defaults.m           - sets the defaults for VBM8
%   cg_vbm8_defs.m               - apply deformations to images
%   cg_vbm8_get_defaults.m       - defaults for VBM8
%   cg_vbm8_longitudinal.m       - VBM8 for longitudinal data
%   cg_vbm8_longitudinal_multi.m - VBM8 for longitudinal data
%   cg_vbm8_longitudinal_multi_run.m - VBM8 for longitudinal data
%   cg_vbm8_run.m                - runtime funtion for VBM8
%   cg_vbm8_tools.m              - wrapper for calling VBM8 utilities
%   cg_vbm8_update.m             - check for new updates
%   cg_vbm8_write.m              - write out VBM8 results
%   spm_vbm8.m                   - toolbox wrapper to call functions
%   tbx_cfg_vbm8.m               - configure VBM8
%
% Utility functions
%   AmapMex.m                    - compilation wrapper for AmapMex.c
%   GBM.m                        - skull-stripping using graph-cut
%   cg_cfg_realign.m             - configuration file for cg_realign.m
%   cg_check_cov.m               - check sample homogeneity across sample
%   cg_cleanup_gwc.m             - use morphological operations to cleanup GM/WM/CSF
%   cg_morph_vol.m               - morphological operations to 3D data
%   cg_realign.m                 - estimation of within modality rigid body movement parameters
%   cg_run_realign_estimate.m    - runtime function for estimate realign
%   cg_run_realign_estwrite.m    - runtime function for estimate and write realign
%   cg_sanlm.m                   - Spatial Adaptive Non Local Means Denoising Filter
%   cg_showslice_all.m           - show 1 slice of all images
%   cg_slice_overlay.m           - wrapper for overlay tool slice_overlay
%   cg_slice_overlay_ui.m        - example for user interface for overlay wrapper cg_slice_overlay.m
%   cg_spmF2x.m                  - transformation of F-maps to P, -log(P), R2 maps
%   cg_spmT2x.m                  - transformation of t-maps to P, -log(P), r or d-maps
%   checkinopt.m                 - check input and options
%   dp.m                         - runtime estimation
%   sanlmMex.m                   - compilation wrapper for sanlmMex.c
%   slice_overlay.m              - overlay tool
%
% Mex- and c-functions
%   Amap.c                       - Adaptive Maximum A Posteriori segmentation
%   Amap.h                       - header for Amap.c
%   AmapMex.c                    - mex-wrapper for Amap 
%   Kmeans.c                     - tree structure k-means algorithm
%   MrfPrior.c                   - estimation of MRF weighting
%   Pve.c                        - partial volume estimaion (PVE)
%   down_cut.c                   - graph-cut functions
%   eikonal3.c                   - eikonal distance calculation for 3D images
%   median3.c                    - median filter for 3D images
%   sanlmMex.c                   - mex-wrapper for sanlm_float.c
%   sanlm_float.c                - Adaptive Non-Local Means Denoising Filter (core functions)
%   vbdist.c                     - voxel-based euclidean distance calculation
%   vollib.c                     - volume convolving functions
%
% Batch functions
%   cg_spm8_batch.m              - batch mode wrapper for spm_jobman for SPM8
%   cg_spm8_batch.sh             - shell script to call matlab batch files from unix
%                                  without gui
%   cg_vbm8_batch.m              - batch mode wrapper for spm_jobman for VBM8
%   cg_vbm8_batch.sh             - shell script to use vbm from unix without gui
%
% Templates/Images
%   Template_?_IXI550_MNI152.nii - Dartel template of 550 subjects from IXI database
%                                  in MNI152 space provided for 6 different iteration steps
% avgT1_Dartel_IXI550_MNI152.nii - average of 550 T1 images of IXI database in MNI152 space
