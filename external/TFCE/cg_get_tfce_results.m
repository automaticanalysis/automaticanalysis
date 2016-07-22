function [SPM,xSPM] = cg_get_tfce_results(varargin)
% Compute a specified and thresholded SPM/PPM following parameter estimation
% FORMAT [SPM,xSPM] = cg_getTFCE;
% Query SPM in interactive mode.
%
% FORMAT [SPM,xSPM] = cg_getTFCE(xSPM);
% Query SPM in batch mode. See below for a description of fields that may
% be present in xSPM input. Values for missing fields will be queried
% interactively.
%
% xSPM      - structure containing SPM, distribution & filtering details
% .swd      - SPM working directory - directory containing current SPM.mat
% .title    - title for comparison (string)
% .Z        - minimum of Statistics {filtered on u and k}
% .n        - conjunction number <= number of contrasts
% .STAT     - distribution {Z, T, X, F or P}
% .df       - degrees of freedom [df{interest}, df{residual}]
% .STATstr  - description string
% .Ic       - indices of contrasts (in SPM.xCon)
% .Im       - indices of masking contrasts (in xCon)
% .pm       - p-value for masking (uncorrected)
% .Ex       - flag for exclusive or inclusive masking
% .u        - height threshold
% .k        - extent threshold {voxels}
% .XYZ      - location of voxels {voxel coords}
% .XYZmm    - location of voxels {mm}
% .S        - search Volume {voxels}
% .R        - search Volume {resels}
% .FWHM     - smoothness {voxels}
% .M        - voxels -> mm matrix
% .iM       - mm -> voxels matrix
% .VOX      - voxel dimensions {mm} - column vector
% .DIM      - image dimensions {voxels} - column vector
% .Vspm     - Mapped statistic image(s)
% .Ps       - uncorrected P values in searched volume (for voxel FDR)
% .Pp       - uncorrected P values of peaks (for peak FDR)
% .Pc       - uncorrected P values of cluster extents (for cluster FDR)
% .uc       - 0.05 critical thresholds for FWEp, FDRp, FWEc, FDRc
% .thresDesc - description of height threshold (string)
%
% Required fields of SPM
%
% xVol   - structure containing details of volume analysed
%
% xX     - Design Matrix structure
%        - (see spm_spm.m for structure)
%
% xCon   - Contrast definitions structure array
%        - (see also spm_FcUtil.m for structure, rules & handling)
% .name  - Contrast name
% .STAT  - Statistic indicator character ('T', 'F' or 'P')
% .c     - Contrast weights (column vector contrasts)
% .X0    - Reduced design matrix data (spans design space under Ho)
%          Stored as coordinates in the orthogonal basis of xX.X from spm_sp
%          (Matrix in SPM99b)  Extract using X0 = spm_FcUtil('X0',...
% .iX0   - Indicates how contrast was specified:
%          If by columns for reduced design matrix then iX0 contains the
%          column indices. Otherwise, it's a string containing the
%          spm_FcUtil 'Set' action: Usually one of {'c','c+','X0'}
% .X1o   - Remaining design space data (X1o is orthogonal to X0)
%          Stored as coordinates in the orthogonal basis of xX.X from spm_sp
%          (Matrix in SPM99b)  Extract using X1o = spm_FcUtil('X1o',...
% .eidf  - Effective interest degrees of freedom (numerator df)
%        - Or effect-size threshold for Posterior probability
% .Vcon  - Name of contrast (for 'T's) or ESS (for 'F's) image
% .Vspm  - Name of SPM image
%
% Evaluated fields in xSPM (input)
%
% xSPM      - structure containing SPM, distribution & filtering details
% .swd      - SPM working directory - directory containing current SPM.mat
% .title    - title for comparison (string)
% .Ic       - indices of contrasts (in SPM.xCon)
% .Im       - indices of masking contrasts (in xCon)
% .pm       - p-value for masking (uncorrected)
% .Ex       - flag for exclusive or inclusive masking
% .u        - height threshold
% .k        - extent threshold {voxels}
% .thresDesc - description of height threshold (string)
%
% In addition, the xCon structure is updated. For newly evaluated
% contrasts, SPM images (spmT_????.{img,hdr}) are written, along with
% contrast (con_????.{img,hdr}) images for SPM{T}'s, or Extra
% Sum-of-Squares images (ess_????.{img,hdr}) for SPM{F}'s.
% 
% The contrast images are the weighted sum of the parameter images,
% where the weights are the contrast weights, and are uniquely
% estimable since contrasts are checked for estimability by the
% contrast manager. These contrast images (for appropriate contrasts)
% are suitable summary images of an effect at this level, and can be
% used as input at a higher level when effecting a random effects
% analysis. (Note that the ess_????.{img,hdr} and
% SPM{T,F}_????.{img,hdr} images are not suitable input for a higher
% level analysis.) See spm_RandFX.man for further details.
%
%__________________________________________________________________________
%
% cg_getTFCE prompts for an SPM and applies thresholds {u & k}
% to a point list of voxel values (specified with their locations {XYZ})
% This allows the SPM be displayed and characterized in terms of regionally 
% significant effects by subsequent routines.
% 
% For general linear model Y = XB + E with data Y, design matrix X,
% parameter vector B, and (independent) errors E, a contrast c'B of the
% parameters (with contrast weights c) is estimated by c'b, where b are
% the parameter estimates given by b=pinv(X)*Y.
% 
% Either single contrasts can be examined or conjunctions of different
% contrasts. Contrasts are estimable linear combinations of the
% parameters, and are specified using the SPM contrast manager
% interface [spm_conman.m]. SPMs are generated for the null hypotheses
% that the contrast is zero (or zero vector in the case of
% F-contrasts). See the help for the contrast manager [spm_conman.m]
% for a further details on contrasts and contrast specification.
% 
% A conjunction assesses the conjoint expression of multiple effects. The
% conjunction SPM is the minimum of the component SPMs defined by the
% multiple contrasts.  Inference on the minimum statistics can be
% performed in different ways.  Inference on the Conjunction Null (one or
% more of the effects null) is accomplished by assessing the minimum as
% if it were a single statistic; one rejects the conjunction null in
% favor of the alternative that k=nc, that the number of active effects k
% is equal to the number of contrasts nc.  No assumptions are needed on
% the dependence between the tests.
%
% Another approach is to make inference on the Global Null (all effects
% null).  Rejecting the Global Null of no (u=0) effects real implies an
% alternative that k>0, that one or more effects are real.   A third
% Intermediate approach, is to use a null hypothesis of no more than u
% effects are real.  Rejecting the intermediate null that k<=u implies an
% alternative that k>u, that more than u of the effects are real.  
% 
% The Global and Intermediate nulls use results for minimum fields which
% require the SPMs to be identically distributed and independent. Thus,
% all component SPMs must be either SPM{t}'s, or SPM{F}'s with the same
% degrees of freedom. Independence is roughly guaranteed for large
% degrees of freedom (and independent data) by ensuring that the
% contrasts are "orthogonal". Note that it is *not* the contrast weight
% vectors per se that are required to be orthogonal, but the subspaces of
% the data space implied by the null hypotheses defined by the contrasts
% (c'pinv(X)). Furthermore, this assumes that the errors are
% i.i.d. (i.e. the estimates are maximum likelihood or Gauss-Markov. This
% is the default in spm_spm).  
%
% To ensure approximate independence of the component SPMs in the case of
% the global or intermediate null, non-orthogonal contrasts are serially
% orthogonalised in the order specified, possibly generating new
% contrasts, such that the second is orthogonal to the first, the third
% to the first two, and so on.  Note that significant inference on the
% global null only allows one to conclude that one or more of the effects
% are real.  Significant inference on the conjunction null allows one to
% conclude that all of the effects are real.
%
% Masking simply eliminates voxels from the current contrast if they
% do not survive an uncorrected p value (based on height) in one or
% more further contrasts.  No account is taken of this masking in the
% statistical inference pertaining to the masked contrast.
% 
% The SPM is subject to thresholding on the basis of height (u) and the
% number of voxels comprising its clusters {k}. The height threshold is
% specified as above in terms of an [un]corrected p value or
% statistic.  Clusters can also be thresholded on the basis of their
% spatial extent. If you want to see all voxels simply enter 0.  In this
% instance the 'set-level' inference can be considered an 'omnibus test'
% based on the number of clusters that obtain.
%
% BAYESIAN INFERENCE AND PPMS - POSTERIOR PROBABILITY MAPS
%
% If conditional estimates are available (and your contrast is a T
% contrast) then you are asked whether the inference should be 'Bayesian'
% or 'classical' (using GRF).  If you choose Bayesian the contrasts are of
% conditional (i.e. MAP) estimators and the inference image is a
% posterior probability map (PPM).  PPMs encode the probability that the
% contrast exceeds a specified threshold.  This threshold is stored in
% the xCon.eidf.  Subsequent plotting and tables will use the conditional
% estimates and associated posterior or conditional probabilities.
% 
% see spm_results_ui.m for further details of the SPM results section.
% see also spm_contrasts.m
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes, Karl Friston & Jean-Baptiste Poline
%
% based on
% $Id: cg_get_tfce_results.m 70 2014-11-13 10:44:52Z gaser $


%-Select SPM.mat & note SPM results directory
%--------------------------------------------------------------------------
if (nargin > 0) 
    xSPM = varargin{1};
end
try
   swd = xSPM.swd;
   sts = 1;
catch
    [spmmatfile, sts] = spm_select(1,'^SPM\.mat$','Select SPM.mat');
    if strcmp(spm('ver'),'SPM12')
        swd = spm_file(spmmatfile,'fpath');
    else
        swd = spm_str_manip(spmmatfile,'H');
    end
end;
if ~sts, SPM = []; xSPM = []; return; end

%-Preliminaries...
%==========================================================================

%-Load SPM.mat
%--------------------------------------------------------------------------
try
    load(fullfile(swd,'SPM.mat'));
catch
    error(['Cannot read ' fullfile(swd,'SPM.mat')]);
end
SPM.swd = swd;


%-Change directory so that relative filenames are valid
%--------------------------------------------------------------------------
cd(SPM.swd);

%-Check the model has been estimated
%--------------------------------------------------------------------------
try
    XYZ  = SPM.xVol.XYZ;
catch
    
    %-Check the model has been estimated
    %----------------------------------------------------------------------
    str = { 'This model has not been estimated.';...
            'Would you like to estimate it now?'};
    if spm_input(str,1,'bd','yes|no',[1,0],1)
         SPM = spm_spm(SPM);
    else
        SPM = []; xSPM = [];
        return
    end
end

xX   = SPM.xX;                      %-Design definition structure
XYZ  = SPM.xVol.XYZ;                %-XYZ coordinates
S    = SPM.xVol.S;                  %-search Volume {voxels}
R    = SPM.xVol.R;                  %-search Volume {resels}
M    = SPM.xVol.M(1:3,1:3);         %-voxels to mm matrix
VOX  = sqrt(diag(M'*M))';           %-voxel dimensions

% check the data and other files have valid filenames
%-----------------------------------------------------------------------
try, SPM.xY.VY     = spm_check_filename(SPM.xY.VY);     end
try, SPM.xVol.VRpv = spm_check_filename(SPM.xVol.VRpv); end
try, SPM.Vbeta     = spm_check_filename(SPM.Vbeta);     end
try, SPM.VResMS    = spm_check_filename(SPM.VResMS);    end
try, SPM.VM        = spm_check_filename(SPM.VM);        end

%-Contrast definitions
%==========================================================================

%-Load contrast definitions (if available)
%--------------------------------------------------------------------------
try
    xCon = SPM.xCon;
catch
    xCon = {};
end

%-Check whether mesh are detected if we use spm12
%--------------------------------------------------------------------------

if exist(fullfile(swd, 'mask.img'))
    file_ext = '.img';
elseif exist(fullfile(swd, 'mask.nii'))
    file_ext = '.nii';
elseif exist(fullfile(swd, 'mask.gii'))
    file_ext = '.gii';
else
    error('No mask file found.');
end

if strcmp(spm('ver'),'SPM12')
    if spm_mesh_detect(SPM.xY.VY)
        mesh_detected = 1;
    else
        mesh_detected = 0;
    end
else
      mesh_detected = 0;
end


%==========================================================================
% - C O N T R A S T S ,   S P M    C O M P U T A T I O N ,    M A S K I N G
%==========================================================================

%-Get contrasts
%--------------------------------------------------------------------------
try
    Ic = xSPM.Ic;
catch
    [Ic,xCon] = spm_conman(SPM,'T',Inf,...
               '    Select contrasts...',' for conjunction',1);
end
if isempty(xCon)
    % figure out whether new contrasts were defined, but not selected
    % do this by comparing length of SPM.xCon to xCon, remember added
    % indices to run spm_contrasts on them as well
    try
        noxCon = numel(SPM.xCon);
    catch
        noxCon = 0;
    end
    IcAdd = (noxCon+1):numel(xCon);
else
    IcAdd = [];
end

nc        = length(Ic);  % Number of contrasts

%-Allow user to extend the null hypothesis for conjunctions
%
% n: conjunction number
% u: Null hyp is k<=u effects real; Alt hyp is k>u effects real 
%    (NB Here u is from Friston et al 2004 paper, not statistic thresh).
%                  u         n      
% Conjunction Null nc-1      1     |    u = nc-n
% Intermediate     1..nc-2   nc-u  |    #effects under null <= u
% Global Null      0         nc    |    #effects under alt  > u,  >= u+1
%----------------------------------+---------------------------------------
if (nc > 1)
  error('No conjunction allowed.');
else
    n = 1;
end

SPM.xCon = xCon;

%-Get contrasts for masking
%--------------------------------------------------------------------------
Im = [];
pm = [];
Ex = [];


%-Create/Get title string for comparison
%--------------------------------------------------------------------------
str  = xCon(Ic).name;
try
    titlestr = xSPM.title;
    if isempty(titlestr)
        titlestr = str;
    end
catch
    titlestr     = spm_input('title for comparison',1,'s',str);
end

titlestr = ['Nonparametric test: ' titlestr];

%-Compute & store contrast parameters, contrast/ESS images, & SPM images
%==========================================================================
SPM.xCon = xCon;
SPM      = spm_contrasts(SPM, unique([Ic, Im, IcAdd]));
xCon     = SPM.xCon;
STAT     = xCon(Ic(1)).STAT;
VspmSv   = cat(1,xCon(Ic).Vspm);

% check that statistic for this contrast was estimated
%--------------------------------------------------------------------------
if ~exist(fullfile(swd, sprintf('T_log_p_%04d%s',Ic,file_ext)))
  strtmp = { 'No TFCE calculation for this contrast found.';...
      'Would you like to estimate it now?'};
  if spm_input(strtmp,1,'bd','yes|no',[1,0],1)
    spm_jobman('interactive','','spm.tools.tfce_estimate');
  end

  SPM = [];
  xSPM = [];
  return
end

%-Degrees of Freedom and STAT string describing marginal distribution
%--------------------------------------------------------------------------
df          = [xCon(Ic(1)).eidf xX.erdf];
try
    stattype = xSPM.stattype;
catch
    stattype = spm_input('Type of statistic','+1','b','TFCE|T',[],1);
end
try
    thresDesc = xSPM.thresDesc;
catch
    thresDesc = spm_input('p value adjustment to control','+1','b','FWE|FDR|none',[],1);
end

switch thresDesc
    case 'FWE' 
        statcorr = 'corrP';
    otherwise
        statcorr = 'P';
end

STAT = stattype;

switch STAT
    case 'TFCE'
        STATstr = sprintf('TFCE_{%.0f}',df(2));
    case 'T'
        STATstr = sprintf('T_{%.0f}',df(2));
    case 'F'
        STATstr = sprintf('F_{%.0f,%.0f}',df(1),df(2));
    case 'P'
        STATstr = sprintf('PPM^{%0.2f}','PPM',df(1));
end

% get # of permutations
try
  n_perm = load(fullfile(swd, sprintf('T_%04d.txt',Ic)));
catch
  n_perm = 0;
end

z_name  = fullfile(swd, sprintf('%s_%04d%s',stattype,Ic,file_ext));
Pz_name = fullfile(swd, sprintf('%s_log_p_%04d%s',stattype,Ic,file_ext));
Pu_name = fullfile(swd, sprintf('%s_log_pFWE_%04d%s',stattype,Ic,file_ext));
Qu_name = fullfile(swd, sprintf('%s_log_pFDR_%04d%s',stattype,Ic,file_ext));

if mesh_detected
    VQu = spm_data_hdr_read(Qu_name);
    VPu = spm_data_hdr_read(Pu_name);
    VPz = spm_data_hdr_read(Pz_name);
    Vz  = spm_data_hdr_read(z_name);
else
    VQu = spm_vol(Qu_name);
    VPu = spm_vol(Pu_name);
    VPz = spm_vol(Pz_name);
    Vz  = spm_vol(z_name);
end

%-Compute SPM
%--------------------------------------------------------------------------
if strcmp(spm('ver'),'SPM12')
    Z   = spm_data_read(Vz,'xyz',XYZ);
    Qu  = spm_data_read(VQu,'xyz',XYZ);
    Pz  = spm_data_read(VPz,'xyz',XYZ);
    Pu  = spm_data_read(VPu,'xyz',XYZ);
else
    Z   = spm_get_data(Vz,XYZ);
    Qu  = spm_get_data(VQu,XYZ);
    Pz  = spm_get_data(VPz,XYZ);
    Pu  = spm_get_data(VPu,XYZ);
end

% convert from -log10
Qu = 1-(10.^-Qu);
Pu = 1-(10.^-Pu);
Pz = 1-(10.^-Pz);

switch thresDesc
    case 'FWE' 
        Zp = Pu;
    case 'FDR' 
        Zp = Qu;
    otherwise
        Zp = Pz;
end

%==========================================================================
% - H E I G H T   &   E X T E N T   T H R E S H O L D S
%==========================================================================

try u = xSPM.u; catch, u = -Inf; end       % height threshold
k   = 0;           % extent threshold {voxels}

%-Use standard FDR 
%--------------------------------------------------------------------------
topoFDR = false;
    
if  mesh_detected
    G = export(gifti(SPM.xVol.G),'patch');
end

%-Height threshold - classical inference
%--------------------------------------------------------------------------
if STAT ~= 'P'
        
    switch thresDesc
        
        case 'FWE' % Family-wise false positive rate
        %------------------------------------------------------------------
        if isinf(u), u = spm_input('p value (FWE)','+0','r',0.05,1,[0,1]); end
        thresDesc = ['p<' num2str(u) ' (' thresDesc ')'];
        u = 1 - u;


        case 'FDR' % False discovery rate
        %------------------------------------------------------------------
        if isinf(u), u = spm_input('p value (FDR)','+0','r',0.05,1,[0,1]); end
        thresDesc = ['p<' num2str(u) ' (' thresDesc ')'];
        u = 1 - u;
        
        case 'none'  % No adjustment
        %------------------------------------------------------------------
        if isinf(u), u = spm_input(['threshold {',STAT,' or p value}'],'+0','r',0.001,1,[0,1]); end
        thresDesc = ['p<' num2str(u) ' (unc.)'];
        u = 1 - u;
        
        otherwise
        %------------------------------------------------------------------
        fprintf('\n');                                                  %-#
        error(sprintf('Unknown control method "%s".',thresDesc));

    end % switch thresDesc
            
end % (if STAT)

%-Calculate height threshold filtering
%--------------------------------------------------------------------------
Q      = find(Zp > u);
if mesh_detected, str = 'vertices'; else str = 'voxels'; end

%-Apply height threshold
%--------------------------------------------------------------------------
Z      = Z(:,Q);
Qu     = Qu(:,Q);
Pz     = Pz(:,Q);
Pu     = Pu(:,Q);
XYZ    = XYZ(:,Q);

if isempty(Q)
    fprintf('\n');                                                      %-#
    sw = warning('off','backtrace');
    warning('SPM:NoVoxels','No %s survive height threshold at u=%0.2g',str,u);
    warning(sw);
end

%-Extent threshold (disallowed for conjunctions)
%--------------------------------------------------------------------------
if ~isempty(XYZ) && nc == 1

    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...extent threshold')  %-#
    
    %-Get extent threshold [default = 0]
    %----------------------------------------------------------------------
    if STAT == 'T'
        try
            k = xSPM.k;
        catch
            k = spm_input(['& extent threshold {' str '}'],'+1','r',0,1,[0,Inf]);
        end
    else
        k = 0;
    end
    
    %-Calculate extent threshold filtering
    %----------------------------------------------------------------------
    if  mesh_detected
        T = false(SPM.xVol.DIM');
        T(XYZ(1,:)) = true;
        A = spm_mesh_clusters(G,T)';
        A = A(XYZ(1,:));
    else
        A = spm_clusters(XYZ);
    end
    Q     = [];
    for i = 1:max(A)
        j = find(A == i);
        if length(j) >= k, Q = [Q j]; end
    end

    % ...eliminate voxels
    %----------------------------------------------------------------------
    Z     = Z(:,Q);
    Qu    = Qu(:,Q);
    Pz    = Pz(:,Q);
    Pu    = Pu(:,Q);
    XYZ   = XYZ(:,Q);
    if isempty(Q)
        fprintf('\n');                                                  %-#
        sw = warning('off','backtrace');
        warning('SPM:NoVoxels','No %s survive extent threshold at k=%0.2g',str,k);
        warning(sw);
    end

else

    k = 0;

end % (if ~isempty(XYZ))

%==========================================================================
% - E N D
%==========================================================================
spm('Pointer','Arrow')

%-Assemble output structures of unfiltered data
%==========================================================================
xSPM   = struct( ...
        'swd',      swd,...
        'title',    titlestr,...
        'Z',        Z,...
        'Qu',       Qu,...
        'Pu',       Pu,...
        'Pz',       Pz,...
        'VQu',      VQu,...
        'VPu',      VPu,...
        'VPz',      VPz,...
        'n',        n,...
        'STAT',     STAT,...
        'df',       df,...
        'STATstr',  STATstr,...
        'Ic',       Ic,...
        'Im',       Im,...
        'pm',       pm,...
        'Ex',       Ex,...
        'u',        u,...
        'k',        k,...
        'XYZ',      XYZ,...
        'XYZmm',    SPM.xVol.M(1:3,:)*[XYZ; ones(1,size(XYZ,2))],...
        'S',        SPM.xVol.S,...
        'R',        SPM.xVol.R,...
        'FWHM',     SPM.xVol.FWHM,...
        'M',        SPM.xVol.M,...
        'iM',       SPM.xVol.iM,...
        'DIM',      SPM.xVol.DIM,...
        'VOX',      VOX,...
        'Vspm',     VspmSv,...
        'n_perm',   n_perm,...
        'thresDesc',thresDesc);

% RESELS per voxel (density) if it exists
%--------------------------------------------------------------------------
try, xSPM.VRpv = SPM.xVol.VRpv; end
try
    xSPM.units = SPM.xVol.units;
catch
    try, xSPM.units = varargin{1}.units; end
end

%-Topology for surface-based inference
%--------------------------------------------------------------------------
if mesh_detected
    xSPM.G     = G;
    xSPM.XYZmm = xSPM.G.vertices(xSPM.XYZ(1,:),:)';
end


