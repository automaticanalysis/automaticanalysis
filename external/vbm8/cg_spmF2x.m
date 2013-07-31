function cg_spmF2x(vargin)
%CG_SPMF2X transformation of F-maps to P, -log(P), R2 maps
%
% The following formulas are used:
%
% --------------------------------
% coefficient of determination R2
% --------------------------------
%
%           F*(n-1)
% R2 = ------------------
%        n-p + F*(n-1)
%
% --------------------------------
% p-value
% --------------------------------
%
% p = 1-spm_Fcdf
%
% --------------------------------
% log p-value
% --------------------------------
%
% -log10(1-P) = -log(1-spm_Fcdf)
%
% For the last case of log transformation this means that a p-value of p=0.99 (0.01)
% is transformed to a value of 2
%
% Examples:
% p-value   -log10(1-P)
% 0.1       1
% 0.05      1.3
% 0.01      2
% 0.001     3
% 0.0001    4
%
% All maps can be thresholded using height and extent thresholds and you can 
% also apply corrections for multiple comparisons based on family-wise error 
% (FWE) or false discovery rate (FDR). You can easily threshold and/or 
% transform a large number of spmF-maps using the same thresholds.
%
% Naming convention of the transformed files:
%   Type_Contrast_Pheight_K
%
%   Type:      P    - p-value
%              logP - log p-value
%              R2   - coefficient of determination
%              F    - F-value
%
%   Contrast:  name used in the contrast manager with replaced none valid 
%              strings
%    
%   Pheight:   p    - uncorrected p-value in % (p<0.05 will coded with "p5")
%              pFWE - p-value with FWE correction in %
%              pFDR - p-value with FDR correction in %
%              
%   K:         extent threshold in voxels
%
%_______________________________________________________________________
% Christian Gaser
% $Id: cg_spmF2x.m 404 2011-04-11 10:03:40Z gaser $

rev = '$Rev: 404 $';

if nargin == 1
    P = char(vargin.data);

    sel = vargin.conversion.sel;

    if isfield(vargin.conversion.threshdesc,'fwe')
        adjustment = 1;
        u0  = vargin.conversion.threshdesc.fwe.thresh;
    elseif isfield(vargin.conversion.threshdesc,'fdr')
        adjustment = 2;
        u0  = vargin.conversion.threshdesc.fdr.thresh;
    elseif isfield(vargin.conversion.threshdesc,'uncorr')
        adjustment = 0;
        u0  = vargin.conversion.threshdesc.uncorr.thresh;
    elseif isfield(vargin.conversion.threshdesc,'none')
        adjustment = -1;
        u0  = -Inf;
    end
    
    if isfield(vargin.conversion.cluster,'fwe')
        extent_FWE = 1;
        pk  = vargin.conversion.cluster.fwe.thresh;
    elseif isfield(vargin.conversion.cluster,'uncorr')
        extent_FWE = 0;
        pk  = vargin.conversion.cluster.uncorr.thresh;
    elseif isfield(vargin.conversion.cluster,'k')
        extent_FWE = 0;
        pk  = vargin.conversion.cluster.k.kthresh;
    else
        extent_FWE = 0;
        pk=0;
    end
end

spm2 = 0;
if strcmp(spm('ver'),'SPM2'), spm2 = 1; end

if nargin < 1
    if spm2
        P = spm_get(Inf,'spmF*.img','Select F-images');
    else
        P = spm_select(Inf,'^spmF.*\.img$','Select F-images');
    end

    sel = spm_input('Convert F value to?',1,'m',...
    '1-p|-log(1-p)|coefficient of determination R^2',1:3, 2);

    %-Get height threshold
    %-------------------------------------------------------------------
    str = 'FWE|FDR|uncorr|none';
    adjustment = spm_input('p value adjustment to control','+1','b',str,[],1);
    
    switch adjustment
    case 1 % family-wise false positive rate
        %---------------------------------------------------------------
        u0  = spm_input('p value (family-wise error)','+0','r',0.05,1,[0,1]);
    case 2 % False discovery rate
        %---------------------------------------------------------------    
        u0  = spm_input('p value (false discovery rate)','+0','r',0.05,1,[0,1]);
    case 0  %-NB: no adjustment
        % p for conjunctions is p of the conjunction SPM
        %---------------------------------------------------------------
        u0  = spm_input(['threshold {T or p value}'],'+0','r',0.001,1);
    otherwise  %-NB: no threshold
        % p for conjunctions is p of the conjunction SPM
        %---------------------------------------------------------------
        u0  = -Inf;
    end

    if adjustment > -1 
        pk     = spm_input('extent threshold {k or p-value}','+1','r',0,1);
    else
        pk = 0;
    end
    if (pk < 1) & (pk > 0)
        extent_FWE = spm_input('p value (extent)','+1','b','uncorrected|FWE corrected',[0 1],1);
    end
end


switch adjustment
case 1 % family-wise false positive rate
    p_height_str = '_pFWE';
case 2 % False discovery rate
    p_height_str = '_pFDR';
case 0 %-NB: no adjustment
    p_height_str = '_p';
otherwise  %-NB: no threshold
    p_height_str = '';
end

for i=1:size(P,1)
    spmF = deblank(P(i,:));
    Vspm = spm_vol(spmF);   
    [pth,nm,xt,vr] = spm_fileparts(spmF);

    SPM_name = fullfile(pth, 'SPM.mat');
    
    % SPM.mat exist?
    if ~exist(SPM_name)
       error('SPM.mat not found')
    end

    if strcmp(nm(1:6),'spmF_0') 
        Ic = str2num(nm(length(nm)-2:length(nm)));
    else
        error('Only spmF_0* files can be used');
    end

    load(SPM_name);
    xCon = SPM.xCon;
    df   = [xCon(Ic).eidf SPM.xX.erdf];
    STAT = xCon(Ic).STAT;
    R    = SPM.xVol.R;          %-search Volume {resels}
    S    = SPM.xVol.S;          %-search Volume {voxels}
    XYZ  = SPM.xVol.XYZ;            %-XYZ coordinates
    FWHM = SPM.xVol.FWHM;
    v2r  = 1/prod(FWHM(~isinf(FWHM)));          %-voxels to resels

    switch adjustment
    case 1 % family-wise false positive rate
    %---------------------------------------------------------------
       u  = spm_uc(u0,df,STAT,R,1,S);

    case 2 % False discovery rate
    %---------------------------------------------------------------
       u  = spm_uc_FDR(u0,df,STAT,1,Vspm,0);

    otherwise  %-NB: no adjustment
    % p for conjunctions is p of the conjunction SPM
    %---------------------------------------------------------------
       if (u0 <= 1) & (u0 > 0)
           u = spm_u(u0,df,STAT);
       else
           u = u0;
       end
    end

    F = spm_get_data(Vspm,XYZ);

    %-Calculate height threshold filtering
    %-------------------------------------------------------------------    
    Q      = find(F > u);

    %-Apply height threshold
    %-------------------------------------------------------------------
    F      = F(:,Q);
    XYZ    = XYZ(:,Q);
    if isempty(Q)
        fprintf('No voxels survive height threshold u=%0.2g\n',u);
    end


    %-Extent threshold
    %-----------------------------------------------------------------------
    if ~isempty(XYZ)
        
       if (pk < 1) & (pk > 0)
            if extent_FWE
                Pk = 1;
                k = 0;
                while (Pk >= pk & k<S)
                    k = k + 1;
                    [Pk Pn] = spm_P(1,k*v2r,u,df,STAT,R,1,S);
                end
                p_extent_str = ['_pkFWE' num2str(pk*100)];
            else
                Pn = 1;
                k = 0;
                while (Pn >= pk & k<S)
                    k = k + 1;
                    [Pk Pn] = spm_P(1,k*v2r,u,df,STAT,R,1,S);
                end
                p_extent_str = ['_pk' num2str(pk*100)];
            end
        else
            k = pk;
            p_extent_str = '';
        end

        %-Calculate extent threshold filtering
        %-------------------------------------------------------------------
        A     = spm_clusters(XYZ);
        Q     = [];
        for i = 1:max(A)
            j = find(A == i);
            if length(j) >= k; Q = [Q j]; end
        end


        % ...eliminate voxels
        %-------------------------------------------------------------------
        F     = F(:,Q);
        XYZ   = XYZ(:,Q);
        if isempty(Q)
            fprintf('No voxels survived extent threshold k=%0.2g\n',k);
        end

    else

        k = 0;

    end % (if ~isempty(XYZ))

    if ~isempty(Q)
    
       switch sel
       case 1
          F2x = 1-spm_Fcdf(F,df);
          F2x_name = 'P_';
       case 2
          F2x = -log10(max(eps,1-spm_Fcdf(F,df)));
          F2x_name = 'logP_';
       case 3
    	  	F2x = (df(2)-1)*F./(df(2) - df(1)+F*(df(2) -1));
		      F2x_name = 'R2_';
       end

       str_num = deblank(xCon(Ic).name);

       % replace spaces with "_" and characters like "<" or ">" with "gt" or "lt"
       str_num(findstr(str_num,' ')) = '_';
       strpos = findstr(str_num,' > ');
       if ~isempty(strpos), str_num = [str_num(1:strpos-1) '_gt_' str_num(strpos+1:end)]; end
       strpos = findstr(str_num,' < ');
       if ~isempty(strpos), str_num = [str_num(1:strpos-1) '_lt_' str_num(strpos+1:end)]; end
       strpos = findstr(str_num,'>');
       if ~isempty(strpos), str_num = [str_num(1:strpos-1) 'gt' str_num(strpos+1:end)]; end
       strpos = findstr(str_num,'<');
       if ~isempty(strpos), str_num = [str_num(1:strpos-1) 'lt' str_num(strpos+1:end)]; end
       str_num = spm_str_manip(str_num,'v');
           
       if u0 > -Inf
           name = [F2x_name str_num p_height_str num2str(u0*100) p_extent_str '_k' num2str(k) '.nii'];
       else
           name = [F2x_name str_num '.nii'];
       end
       fprintf('Save %s\n', name);
    
       out = deblank(fullfile(pth,name));

       %-Reconstruct (filtered) image from XYZ & F pointlist
       %-----------------------------------------------------------------------
       Y      = zeros(Vspm.dim(1:3));
       OFF    = XYZ(1,:) + Vspm.dim(1)*(XYZ(2,:)-1 + Vspm.dim(2)*(XYZ(3,:)-1));
       Y(OFF) = F2x;

       VO = Vspm;
       VO.fname = out;
       if spm2
            VO.dim(4) = spm_type('float32');
       else
            VO.dt = [spm_type('float32') spm_platform('bigend')];
       end
       spm_write_vol(VO,Y);
    
    end
end
