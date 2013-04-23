function cg_spmT2x(vargin)
%CG_SPMT2X transformation of t-maps to P, -log(P), r or d-maps
%
% The following formulas are used:
%
% --------------------------------
% correlation coefficient:
% --------------------------------
%          sign(t)
% r = ------------------
%            df
%     sqrt(------ + 1)
%           t*t
%
% --------------------------------
% effect-size
% --------------------------------
%            2t
% d = ----------------
%         sqrt(df)
%
% --------------------------------
% p-value
% --------------------------------
%
% p = 1-spm_Tcdf
%
% --------------------------------
% log p-value
% --------------------------------
%
% -log10(1-P) = -log(1-spm_Tcdf)
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
% transform a large number of spmT-maps using the same thresholds.
%
% Naming convention of the transformed files:
%   Type_Contrast_Pheight_Pextent_K_Neg
%
%   Type:      P    - p-value
%              logP - log p-value
%              R    - correlation coefficient
%              D    - effect size
%              T    - t-value
%
%   Contrast:  name used in the contrast manager with replaced none valid 
%              strings
%    
%   Pheight:   p    - uncorrected p-value in % (p<0.05 will coded with "p5")
%              pFWE - p-value with FWE correction in %
%              pFDR - p-value with FDR correction in %
%              
%   Pextent:   pk    - uncorr. extent p-value in % (p<0.05 coded with "p5")
%              pkFWE - extent p-value with FWE correction in %
%
%   K:         extent threshold in voxels
%
%   Neg:       image also shows thresholded inverse effects (e.g. neg. 
%              values) 
%_______________________________________________________________________
% Christian Gaser
% $Id: cg_spmT2x.m 428 2011-09-29 05:59:24Z gaser $

rev = '$Rev: 428 $';

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
    
    if isfield(vargin.conversion.cluster,'fwe2')
        extent_FWE = 1;
        pk  = vargin.conversion.cluster.fwe2.thresh;
        noniso = vargin.conversion.cluster.fwe2.noniso;
    elseif isfield(vargin.conversion.cluster,'uncorr2')
        extent_FWE = 0;
        pk  = vargin.conversion.cluster.uncorr2.thresh;
        noniso = vargin.conversion.cluster.uncorr2.noniso;
    elseif isfield(vargin.conversion.cluster,'k')
        extent_FWE = 0;
        pk  = vargin.conversion.cluster.k.kthresh;
        noniso = vargin.conversion.cluster.k.noniso;
    elseif isfield(vargin.conversion.cluster,'En')
        extent_FWE = 0;
        pk = -1;
        noniso = vargin.conversion.cluster.En.noniso;
    else
        extent_FWE = 0;
        pk = 0;
        noniso = 0;
    end

    neg_results = vargin.conversion.inverse;
    
end

spm2 = 0;
if strcmp(spm('ver'),'SPM2'), spm2 = 1; end

if nargin < 1
    if spm2
        P = spm_get(Inf,'spmT*.img','Select T-images');
    else
        P = spm_select(Inf,'^spmT.*\.img$','Select T-images');
    end
    sel = spm_input('Convert t value to?',1,'m',...
    '1-p|-log(1-p)|correlation coefficient cc|effect size d|apply thresholds without conversion',1:5, 2);

    %-Get height threshold
    %-------------------------------------------------------------------
    str = 'FWE|FDR|uncorr|none';
    adjustment = spm_input('p value adjustment to control','+1','b',str,[1 2 0 -1],1);
    
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
        pk = spm_input('extent threshold {k or p-value}','+1','r',0,1);
    else
        pk = 0;
    end
    if (pk < 1) & (pk > 0)
        extent_FWE = spm_input('p value (extent)','+1','b','uncorrected|FWE corrected',[0 1],1);
    end

    neg_results = spm_input('Show also inverse effects (e.g. neg. values)','+1','b','yes|no',[1 0],2);

    if pk ~= 0
        noniso = spm_input('Correct for non-isotropic smoothness?','+1','b','no|yes',[0 1],2);
    else
        noniso = 0;
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
    spmT = deblank(P(i,:));
    Vspm = spm_vol(spmT);   
    [pth,nm,xt,vr] = spm_fileparts(spmT);

    SPM_name = fullfile(pth, 'SPM.mat');
    
    % SPM.mat exist?
    if ~exist(SPM_name)
       error('SPM.mat not found')
    end

    if strcmp(nm(1:6),'spmT_0') 
        Ic = str2num(nm(length(nm)-2:length(nm)));
    else
        error('Only spmT_0* files can be used');
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

    if ~isfield(SPM.xVol,'VRpv')
        noniso = 0;
    end

    if noniso
        SPM.xVol.VRpv = spm_vol(fullfile(pth,SPM.xVol.VRpv.fname));
    end

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

    Z = spm_get_data(Vspm,XYZ);

    %-Calculate height threshold filtering
    %-------------------------------------------------------------------    
    if neg_results
        Q      = find((Z > u) | (Z < -u));
    else
        Q      = find(Z > u);
    end

    %-Apply height threshold
    %-------------------------------------------------------------------
    Z      = Z(:,Q);
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
        elseif (pk < 0)
            k = 0;
            [P2 Pn2 Em En] = spm_P(1,k,u,df,STAT,R,1,S);
            k = ceil(En/v2r);
            p_extent_str = '_En';
        else
            k = pk;
            p_extent_str = '';
        end
        
        %-Calculate extent threshold filtering
        %-------------------------------------------------------------------
        A     = spm_clusters(XYZ);
        if noniso
            fprintf('Use local RPV values to correct for non-stationary of smoothness.\n');

            Q     = [];
            [N2 Z2 XYZ2 A2 L2] = spm_max(Z,XYZ);
            if max(A) ~= max(A2)
                 error('Number of clusters does not much in spm_max and spm_clusters!');
            end

            for i2 = 1:max(A2)
                %-Get LKC for voxels in i-th region
                %----------------------------------------------------------
                LKC  = spm_get_data(SPM.xVol.VRpv,L2{i2});
                
                %-Compute average of valid LKC measures for i-th region
                %----------------------------------------------------------
                valid = ~isnan(LKC);
                if any(valid)
                    LKC = sum(LKC(valid)) / sum(valid);
                else
                    LKC = v2r; % fall back to whole-brain resel density
                end
                
                %-Intrinsic volume (with surface correction)
                %----------------------------------------------------------
                IV   = spm_resels([1 1 1],L2{i2},'V');
                IV   = IV*[1/2 2/3 2/3 1]';
                k_noniso = IV*LKC/v2r;

                % find corresponding cluster in spm_clusters if cluster exceeds threshold
                if k_noniso >= k
                    ind2 = find(A2==i2);
                    for i = 1:max(A)
                        j = find(A == i);
                        if length(j)==N2(ind2)
                            if any(ismember(XYZ2(:,ind2)',XYZ(:,j)','rows'))
                                Q = [Q j];
                                break
                            end
                        end
                    end
                end
                
            end
        else
            Q     = [];
            for i = 1:max(A)
                j = find(A == i);
                if length(j) >= k; Q = [Q j]; end
            end
        end


        % ...eliminate voxels
        %-------------------------------------------------------------------
        Z     = Z(:,Q);
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
          t2x = 1-spm_Tcdf(Z,df(2));
          t2x_name = 'P_';
       case 2
          t2x = -log10(max(eps,1-spm_Tcdf(Z,df(2))));
          % find neg. T-values
          ind_neg = find(Z<0);
          if ~isempty(ind_neg)
              t2x(ind_neg) = log10(max(eps,spm_Tcdf(Z(ind_neg),df(2))));
          end
          t2x_name = 'logP_';
       case 3
          t2x = sign(Z).*(1./((df(2)./((Z.*Z)+eps))+1)).^0.5;
          t2x_name = 'R_';
       case 4
          t2x = 2*Z/sqrt(df(2));
          t2x_name = 'D_';
       case 5
          t2x = Z;
          t2x_name = 'T_';
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
    
       if neg_results
            neg_str = '_bi'; 
       else
            neg_str = '';
       end
       
       if u0 > -Inf
           name = [t2x_name str_num p_height_str num2str(u0*100) p_extent_str '_k' num2str(k) neg_str '.nii'];
       else
           name = [t2x_name str_num '.nii'];
       end
       fprintf('Save %s\n', name);
    
       out = deblank(fullfile(pth,name));

       %-Reconstruct (filtered) image from XYZ & Z pointlist
       %-----------------------------------------------------------------------
       Y      = zeros(Vspm.dim(1:3));
       OFF    = XYZ(1,:) + Vspm.dim(1)*(XYZ(2,:)-1 + Vspm.dim(2)*(XYZ(3,:)-1));
       Y(OFF) = t2x;

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
