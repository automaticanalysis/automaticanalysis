% TFCE Toolbox
% Version  67  (TFCE)  2014-04-15
% __________________________________________________________________________
% Copyright (C) 2013 Christian Gaser christian.gaser@uni-jena.de
%
% $Id: Contents_info.txt 62 2013-10-01 14:29:17Z gaser $
% ==========================================================================
% Description
% ==========================================================================
% This toolbox is a an extensions to SPM5/SPM8 (Wellcome Department of Cognitive 
% Neurology) to provide non-parametric statistics based on threshold-free
% cluster enhancement (TFCE). It is developed by Christian Gaser (University of 
% Jena, Department of Psychiatry) and is available to the scientific 
% community under the terms of the GNU General Public License.
%
% General files
%   spm_TFCE.m                  - GUI
%   TFCE.man                    - notes on TFCE toolbox
%   CHANGES.txt                 - changes in revisions
%
% TFCE functions
%   cg_tfce_estimate.m          - TFCE core function
%   cg_tfce_list.m              - list result table
%   cg_tfce_results.m           - call results
%   cg_get_tfce_results.m       - get TFCE information
%   cg_tfce_update.m            - get TFCE update
%   cg_tfce_progress.m          - display progress and remaining time
%   tbx_cfg_tfce_estimate.m     - toolbox function
%   snpm_P_FDR.m                - FDR correction from SnPMs
%
%
% Mex- and c-functions
%   cg_glm_get_Beta_ResSS.c     - GLM estimation
%   tfceMex.c                   - TFCE-transformation 
%   spm_*.c                     - spm functions for volume mapping 
%