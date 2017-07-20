function [P] = snpm_P_FDR(Z,df,STAT,n,Ps)
% Returns the corrected FDR P value 
% FORMAT [P] = snpm_P_FDR(Z,df,STAT,n,Ps)
%
% Z     - height (minimum of n statistics)
%         If empty, find all FDR P values for Ps vector (in which case
%         STAT is ignored, as Ps must be P-values already).
% df    - [df{interest} df{error}]
% STAT  - Statisical field
%		'Z' - Gaussian field
%		'T' - T - field
%		'X' - Chi squared field
%		'F' - F - field
%		'P' - P-values
% n     - number of component SPMs in conjunction
% Ps    - Vector of sorted (ascending) p-values in search volume
%
% P     - corrected FDR   P value
%
%___________________________________________________________________________
%
% The Benjamini & Hochberch (1995) False Discovery Rate (FDR) procedure
% finds a threshold u such that the expected FDR is at most q.  snpm_P_FDR
% returns the smallest q such that Z>u. 
%
% If abs(n) > 1, a P-value for a minimum of n values of the statistic 
% is returned.  If n>0, the P-value assesses the conjunction null
% hypothesis of one or more of the null hypotheses being true.  If n<0,
% then the P-value assess the global null of all nulls being true. 
%
% FDR Background
%
% For a given threshold on a statistic image, the False Discovery Rate
% is the proportion of suprathreshold voxels which are false positives.
% Recall that the thresholding of each voxel consists of a hypothesis
% test, where the null hypothesis is rejected if the statistic is larger
% than threshold.  In this terminology, the FDR is the proportion of
% rejected tests where the null hypothesis is actually true.
%
% A FDR proceedure produces a threshold that controls the expected FDR
% at or below q.  The FDR adjusted p-value for a voxel is the smallest q
% such that the voxel would be suprathreshold.
%
% In comparison, a traditional multiple comparisons proceedure
% (e.g. Bonferroni or random field methods) controls Familywise Error
% rate (FWER) at or below alpha.  FWER is the *chance* of one or more
% false positives anywhere (whereas FDR is a *proportion* of false
% positives).  A FWER adjusted p-value for a voxel is the smallest alpha
% such that the voxel would be suprathreshold. 
%
% If there is truely no signal in the image anywhere, then a FDR
% proceedure controls FWER, just as Bonferroni and random field methods
% do. (Precisely, controlling E(FDR) yields weak control of FWE).  If
% there is some signal in the image, a FDR method should be more powerful
% than a traditional method.
%
% For careful definition of FDR-adjusted p-values (and distinction between
% corrected and adjusted p-values) see Yekutieli & Benjamini (1999).
%
%
% References
%
% Benjamini & Hochberg (1995), "Controlling the False Discovery Rate: A
% Practical and Powerful Approach to Multiple Testing". J Royal Stat Soc,
% Ser B.  57:289-300.
%
% Benjamini & Yekutieli (2001), "The Control of the false discovery rate
% in multiple testing under dependency". To appear, Annals of Statistics.
% Available at http://www.math.tau.ac.il/~benja 
%
% Yekutieli & Benjamini (1999). "Resampling-based false discovery rate
% controlling multiple test procedures for correlated test
% statistics".  J of Statistical Planning and Inference, 82:171-196.
%___________________________________________________________________________
% Based on FIL spm_P_FDR.m     2.6 Thomas Nichols 04/07/02
% @(#)snpm_P_FDR.m	3.2 05/02/23
%	$Id: snpm_P_FDR.m 50 2011-06-21 15:01:54Z gaser $	


if n>0
  Cnj_n = 1;       % Inf on Conjunction Null
else
  Cnj_n = abs(n);  % Inf on Global Null
end

if isempty(Z), 
  AllP = 1;
  Z = Ps;
  STAT = 'P';
else
  AllP = 0;
end

% Set Benjamini & Yeuketeli cV for independence/PosRegDep case
%-----------------------------------------------------------------------
cV   = 1; 
S    = length(Ps);

% Calculate p value of Z
%-----------------------------------------------------------------------
if      STAT == 'Z'
  PZ = (1 - spm_Ncdf(Z)).^Cnj_n;
elseif  STAT == 'T'
  PZ = (1 - spm_Tcdf(Z,df(2))).^Cnj_n;
elseif  STAT == 'X'
  PZ = (1 - spm_Xcdf(Z,df(2))).^Cnj_n;
elseif  STAT == 'F'
  PZ = (1-spm_Fcdf(Z,df)).^Cnj_n;
elseif  STAT == 'P'
  PZ = Z;
end

%-Calculate FDR p values
%-----------------------------------------------------------------------
% If Z is a value in the statistic image, then the adjusted p-value
% defined in Yekutieli & Benjamini (1999) (eqn 3) is obtained.  If Z
% isn't a value in the image, then the adjusted p-value for the next
% smallest statistic value (next largest uncorrected p) is returned.

%-"Corrected" p-values
%-----------------------------------------------------------------------
Qs    = Ps*S./(1:S)'*cV;    

% -"Adjusted" p-values
%-----------------------------------------------------------------------
P     = zeros(size(Z));
if AllP
  P(S) = Qs(S);
  for i = S-1:-1:1
    P(i) = min(Qs(i),P(i+1));
  end
else
  for i = 1:length(Z)
    
    % Find PZ(i) in Ps, or smallest Ps(j) s.t. Ps(j) >= PZ(i)
    %---------------------------------------------------------------
    I = min(find(Ps>=PZ(i)));  
    
    % -"Adjusted" p-values
    %---------------------------------------------------------------
    if isempty(I)
      P(i) = 1;
    else
      P(i) = min(Qs(I:S));
    end
  end
end
