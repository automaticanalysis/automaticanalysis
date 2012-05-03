function I = GLM_Flex(I,DD)
%%% This is the main analysis script.  
%%% Go to: http://nmr.mgh.harvard.edu/harvardagingbrain/People/AaronSchultz/Aarons_Scripts.html
%%% for more information on this script and how to use it.
%%%
%%% Inputs: This is the I structure that is passed to GLM_Flex to run the
%%% analyses.
%%%
%%% I.OutputDir = pwd;
%%% I.F = [];
%%% I.Scans = [];
%%% I.RemoveOutliers = 0;
%%% I.DoOnlyAll = 0;
%%% I.CompOpt = 0;
%%% I.minN = 2;
%%% I.minRat = 0;
%%% I.Thresh = [];
%%%
%%% Written by Aaron Schultz (aschultz@martinos.org), in collaboration with
%%% Donald McLaren.
%%%
%%% Copyright (C) 2011,  Aaron P. Schultz
%%%
%%% Supported in part by the NIH funded Harvard Aging Brain Study (P01AG036694) and NIH R01-AG027435 
%%%
%%% This program is free software: you can redistribute it and/or modify
%%% it under the terms of the GNU General Public License as published by
%%% the Free Software Foundation, either version 3 of the License, or
%%% any later version.
%%% 
%%% This program is distributed in the hope that it will be useful,
%%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%% GNU General Public License for more details.

if ~isfield(I,'DoOnlyAll')
    I.DoOnlyAll = 0;
end

if ~isfield(I,'CompOpt')
    I.CompOpt = 0;
end

%%% Setup Defaults if not Specified
if ~isfield(I,'OutputDir')
   I.OutputDir = pwd;
else
    if exist(I.OutputDir)==0
        mkdir(I.OutputDir);
    end
    cd(I.OutputDir);
end

if ~isfield(I,'minN')
    I.minN = 2;
end
if ~isfield(I,'minRAT')
    I.minRAT = 0;
end

if ~isfield(I,'RemoveOutliers');
   I.RemoveOutliers = 0;
end

if ~isfield(I,'writeIni');
   I.writeIni = 0;
end

if ~isfield(I,'writeFin');
   I.writeFin = 1;
end

if ~isfield(I,'writeoo');
   I.writeoo = 1;
end

if ~isfield(I,'writeI');
   I.writeI = 1;
end

if ~isfield(I,'KeepResiduals');
   I.KeepResiduals = 0;
end

if ~isfield(I,'writeRes');
   I.writeRes = 1;
end

if nargin==1
    ok = 0;
    if exist([I.OutputDir filesep 'IniDat.mat']);
        load([I.OutputDir filesep 'IniDat.mat'])
        try
            tmp = char(scns)==char(I.Scans);
            if mean(tmp(:)) == 1
                [trash v] = openIMG(I.Scans{1});
                if ~isfield(v,'dim');
                    v.dim = size(trash);
                end
                clear trash;
                I.v = v;
                ok = 1;
            else
                ok=0;
            end
        catch
            ok=0;
        end
    end
    if ok == 0;
        if ~isfield(I,'Scans');
            I.Scans = [];
            for ii = 1:length(I.F.IN.N_subs);
                if isempty(I.F.IN.Within)
                    n = I.F.IN.N_subs(ii);
                    g = spm_select(n,'img',['Select Images for Group #' num2str(ii)]);
                    I.Scans = [I.Scans; g];
                else
                    error('this option is not configured');
                end
                
            end
        end
        
        if isfield(I,'Scans');
            fprintf('Reading In Image Files:\n\n');
            s2 = sum(I.F.IN.N_subs).*prod(I.F.IN.Within);          
            [trash v] = openIMG(I.Scans{1});
            v.dim = size(trash);
            I.v = v;
            clear trash
            
            [DD VY] = openIMG(char(I.Scans));
            DD = double(DD);
            
            ss = size(DD);
            if numel(ss)>2
                DD = reshape(DD, prod(ss(1:3)),ss(4))';
            end
            if isfield(I,'NoMiss');
                if I.NoMiss == 1;
                    tmp = mean(diff(DD));
                    DD(:,(tmp==0 | isnan(tmp))) = NaN;
                end
            end
            scns = I.Scans;
            if I.writeIni == 1;
                save IniDat.mat DD scns -v7.3;
            end
            
            fprintf('\n\n');
        end
    end
else
    %DD was passed in as a preconfigured matrix.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(I,'Transform')
    if isfield(I.Transform,'FisherZ')
        if I.Transform.FisherZ == 1;
            DD = atanh(DD);
            disp('FisherZ transform worked!');
        end
    end
    if isfield(I.Transform,'AdjustR2')
        if I.Transform.AdjustR2 == 1;
            DD = atanh(sqrt(DD));
            disp('R2 adjust worked!');
        end
    end
    if isfield(I.Transform,'ZScore')
        if I.Transform.ZScore == 1;
            DD = zscore(DD,0,1);
            disp('Normalization worked!');
        end
    end
end

if isfield(I,'Mask');
    mh = spm_vol(I.Mask);
    msk = resizeVol(mh,I.v);
    DD(:,(msk(:)==0))=NaN;
end

%%% Get the Design Matrix
X1 = I.F.XX;
X = X1; X(isnan(X))=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Setup Objects for Output Data
ss = I.v.dim;
NN = nan(ss);
ResMS = cell(1,numel(I.F.err));
ResMS(:) = {NN};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MasterRes = nan(size(DD));
%%% Get the original full list of voxel indices;
OrigIndex = 1:size(DD,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save I.mat I -v7.3;

%%% Initialize a counter for the subsequent for loops
BigCount = 1;
%%% Screen Missing Data for Within Subject Factors;
if ~isempty(I.F.IN.Within);
   Subs = I.F.FF{size(I.F.FF,1),size(I.F.FF,1)};
   Subs(isnan(Subs))=0;
    for ii = 1:size(Subs,2);
        i2 = find(Subs(:,ii)==1);
        DD(i2,isnan(sum(DD(i2,:)))) = NaN;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
if(size(X,1)~=size(DD,1))
   error('Design Matrix does not match input volumes'); 
end

DD(:,find(sum(DD)==0))=NaN;
%%% Get the number of observations for each Voxel
counts = sum(~isnan(DD));
SS = size(DD,1);
oCy = zeros(SS,SS);
oN  = zeros(SS,SS);
First = 0;
%%% Begin looping through the possible N-sizes 
persisText
for ll = size(X,1):-1:I.minN
    for z = 1
    if I.CompOpt==1
        if exist([pwd filesep 'Comp'])==0
            mkdir Comp;
            cd Comp
        else
            delete Comp/*;
            cd Comp
        end
        
        F = I.F;
        F.Scans = I.Scans;
        RunSPM8ana(F);
        
        I.CompOpt = 0;

        I = GLM_Flex(I);
        return
    end
    end
    %%% Find the indices where the number of observations matches the
    %%% current loop.
    CurrIndex = find(counts==ll);

    if isempty(CurrIndex)
        continue
    end
      
    %%% This Section computes different combinations of N size ll across
    %%% the groups specified in the design matrix.
    if I.DoOnlyAll ~= 1
        rr = ~isnan(DD(:,CurrIndex)).*repmat((1:size(X,1))',1,size(DD(:,CurrIndex),2));
        
        rr = rr';
        [rr II] = sortrows(rr,(1:size(rr,2))*-1);
        
        
        tt = sum(abs(diff(rr~=0)),2)>0;
        ind = find(tt==1);
        ch = [[1; ind+1] [ind; numel(II)]];
        ch = ch((diff(ch,1,2)>=0),:);
        uni = rr(ch(:,1),:);
        
        CurrIndex = CurrIndex(II);
    else
        ch = [1 numel(CurrIndex)];
        uni = 1:size(DD,1);
    end
    
    persisText(['Processing Set #' num2str(ll)]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Loop through specific designs for ll observations
    %disp(numel(ch));
    for ii = 1:size(ch,1)
        %%% Get the observation index for the current design
        Xind = uni(ii,:);
        Xind = Xind(Xind>0);
        
        %%% Subset the design matrix
        xx = X(Xind,:);
        if BigCount==1
            I.X = xx;
        end
        
        %%% Make sure there is enough data across conditions.
        %keyboard
        ttxx = xx(:,setdiff(1:numel(I.F.name), contains('^S.*',I.F.name)));
        ttxx = ttxx(:, ~I.F.CovarCols(setdiff(1:numel(I.F.name), contains('^S.*',I.F.name))));
        if mean(sum(ttxx)>=I.minN)~=1
            continue;
        end
        if  min(sum(ttxx))/max(sum(ttxx)) < I.minRAT;
            continue;
        end
        
        
        %%% get the global index of the voxels being analyzed
        vec = CurrIndex((ch(ii,1):ch(ii,2)));
        %%% Get the correponding data to be analyzed
        Y = DD(Xind,vec);
        
        %%% Estimate the DFs for the full model
        if size(xx,2)==1;
            df1 = 1;
        else
            [z1 z2 z3] = svd(xx);
            tol = max(size(xx))*max(abs(diag(z2)))*eps;
            df1 = sum(diag(z2)>tol);
        end
        df2 = size(xx,1)-df1;
        
        %%% Run a basic GLM
        pv = pinv(xx); 
        pv(abs(pv)<eps*(size(xx,1)))=0;
        beta  = pv*Y;                 
        pred = xx*beta;
        res   = Y-pred;                     
        ResSS = sum(res.^2);                    
        MSE = ResSS./df2;

        %%% Outlier Detection Compute Cook's D
        if I.RemoveOutliers == 1;

            e2 = res.^2;
            p = size(xx,2);
            bc = pinv(xx'*xx); 
            bc(abs(bc)<eps*(size(xx,1)))=0; 
            tmp = diag(xx*bc*xx');
            CD = (e2./(p*repmat(MSE,size(e2,1),1))) .* (repmat(tmp./((1-tmp).^2),1,size(e2,2)));

            if isfield(I,'Thresh') && ~isempty(I.Thresh)
                thresh = icdf('F',I.Thresh,df1,size(xx,1)-df1-1);
            else
                thresh = icdf('F',.5,df1,size(xx,1)-df1-1);
                %disp(thresh);
            end
            
            [mm,ith] = max(CD);
            a = find(mm>thresh);
            
            if ~isempty(a)
                
                scans = Xind(ith(a));
                cols = vec(a);
                vecind = sub2ind(size(DD),scans', cols');
                DD(vecind) = NaN;
                counts(cols) = ll-1;
                
                if ~isempty(I.F.IN.Within);
                    Subs = I.F.FF{size(I.F.FF,1),size(I.F.FF,1)};
                    Subs(isnan(Subs))=0;
                    for jj = 1:size(Subs,2);
                        i2 = find(Subs(:,jj)==1);
                        i3 = find(isnan(sum(DD(i2,:))));
                        DD(i2,i3) = NaN;
                    end
                    counts(cols)=ll-(prod(I.F.IN.Within)-1);
                end
            end
            a = find(mm<thresh);
        else
            a = 1:numel(vec);
        end
        
        %%% If there are no voxels free of outliers, continue
        if isempty(a)
            continue;
        end
        
        First = First+1;
        
        if First
            h = I.v;
            h.fname = 'AllMask.nii';
            tmp = nan(h.dim);
            tmp(vec(a))=1;
            h.dt = [2 0];
            spm_write_vol(h,tmp);
        end
        
        %%% Subset data to only those voxels without outliers
        Y = Y(:,a);
        
        if numel(I.F.Vi)>1
            %%% Compute significance of full model for all voxels
            PSS = sum(pred(:,a).^2);
            ESS = sum(res(:,a).^2);
            FF = (PSS./df1) ./ (ESS./df2);
            UF = icdf('f',1-.001,df1,df2);

            %%% Get the index of voxels where the full model p is less than 0.001
            mv = find(FF>UF);
            cn = numel(mv);
            
            %%% Pool variance across voxels
            if ~isempty(mv)
                q  = spdiags(sqrt(df2./ResSS(a(mv))'),0,cn,cn);
                YY = Y(:,mv)*q;
                Cy = (YY*YY');
                
                oCy(Xind,Xind) = oCy(Xind,Xind)+Cy;
                oN(Xind,Xind) = oN(Xind,Xind)+cn;
            else
                Cy = [];
            end
            Cy = oCy(Xind,Xind) ./ oN(Xind,Xind) ;

            %%% Subset the variance/covariance partitions for the current model
            Vi = I.F.Vi;
            for jj = 1:length(Vi)
                Vi{jj} = I.F.Vi{jj}(Xind,Xind);
            end
            
            %%% Estimate population Var/Covar using ReML estimation.
            try
                xxW = xx(:,1:length(I.F.name));
                [V h] = spm_reml(Cy,xxW,Vi);
            catch
                disp('No go');
                continue;
            end
            V = V*size(xx,1)/trace(V);
        else
            V = eye(size(xx,1));
        end
        
        %%% Store the data/results from before Variance Correction.
        oo.Xind{BigCount} = Xind;
        oo.vec{BigCount} = vec(a);
        oo.Y{BigCount} = Y;
        oo.X{BigCount} = xx;
        %oo.beta1{BigCount} = beta(:,a); %%% Just changed this - Sept. 8, 2011
        
        %%% Create the W mixing matrix to correct for unequal variance and
        %%% dependence between observations
        if numel(I.F.Vi)==1
            W = V;
        else
            W     = full(spm_sqrtm(spm_inv(V)));
            W     = W.*(abs(W) > 1e-6);
        end        
        %%% Store these
        oo.V{BigCount} = V;
        oo.W{BigCount} = W;
        
        %%% Create the alternate design matrix
        aX = W*xx;
        if BigCount == 1;
           I.X = xx;
           I.aX = aX;
        end
        %%% Create the Corrected Data
        Yp = W*Y;
        
        pv = pinv(aX); 
        pv(abs(pv)<eps*(size(xx,1)))=0;
        betap = pv*Yp;
        
        %%% Just changed this - Sept. 8, 2011
        oo.beta1{BigCount} = betap;
        %%% Estimate the partitioned error variance for the current model
        F = I.F;
        F.WX = aX;
        F.XX = I.F.XX(Xind,:);
        F.ind = Xind;
        
        try
            if First && I.writeRes==1
                %keyboard;
                F.writeRes = 1;
                F.v = v;
                F.vec = vec(a);
                F = estimateError(Yp,F);   
                %I.eDF = F.ErrorDF;
                I.FWHM = F.FWHM;
                %I.ReslInfo = F.ReslInfo;
                if ~I.KeepResiduals
                    delete ResAll*.nii ResAll*.mat;
                end
            else
                F = estimateError(Yp,F);
            end
        catch
            warning('The error term was not estimated for this data set');
            continue;
        end
        %%% Store it
        oo.ESS{BigCount} = F.ErrorSS;
        oo.EDF{BigCount} = F.ErrorDF;
        oo.eCons{BigCount} = F.eCons;

        %%% Put the data into images and data structures
        for jj = 1:length(F.ErrorSS);
            ResMS{jj}(vec(a)) = F.ErrorSS{jj}./F.ErrorDF(jj);
        end
        
        NN(vec(a)) = numel(Xind);
        BigCount = BigCount+1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    if I.DoOnlyAll == 1        
        break;
    end
end

if ~exist('oo','var')
   error('Something went wrong.  No voxels were anlyzed.  This is most commonly cause by a mis-specified I.minN or a problem with performing the variance/covariance correction');
end

fprintf('\n\nWriting out mat files and image files...');
if I.writeFin == 1;
    save FinDat.mat DD -v7.3;
end

v = I.v;
v.dt = [64 0];
% v.descrip = [];

%%%
for ll = 1:length(ResMS);
    tm = 'ResMS_00';
    n = num2str(ll);
    tm((end-length(n))+1:end) = n;
    writeIMG(v,ResMS{ll},[tm '.nii']);
end

writeIMG(v,NN,'NN.nii');
writeIMG(v,NN>0,'mask.nii');

if isfield(I,'Cons')  
   [I oo] = GLM_Flex_Contrasts(I,1:length(I.Cons),oo); 
end

if I.writeI == 1;
    save I.mat I -v7.3;
end

if I.writeoo == 1;
    save oo.mat oo -v7.3;
end

