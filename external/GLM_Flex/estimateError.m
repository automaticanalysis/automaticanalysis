function [F RES] = estimateError(y,F)
%%% This script will compute the error terms specified in the Design
%%% Structure F.  
%%%
%%% Inputs:
%%% y = a column vector or matrix of data mated to the design specified in F
%%%
%%% F = Design structure from CreateDesign.m
%%%.
%%%
%%% Outputs:
%%% The computed statistics are put in the data structure F and returned.
%%%
%%%
%%% Written by Aaron Schultz (aschultz@martinos.org) 
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

if any(F.CovarCols) && ~isempty(F.IN.Within)
    error('You are trying to run a repeated measures ANCOVA.  This analysis is not supported at this time.');
end

if ~isfield(F,'writeRes');
    F.writeRes=0;
end

if isfield(F,'ind');
    Xind = F.ind;
else
    Xind = 1:size(F.XX,1);
end

if size(F.XX,1) == numel(Xind)
	xx = F.XX;
else
	xx = F.XX(Xind,:);
end

if isfield(F,'WX');
	if size(F.WX,1) == numel(Xind)
    	x = F.WX;
    else
    	x = F.WX(Xind,:);
	end
else
    x = xx;
    x(isnan(xx))=0;
end

wh = find(sum(x)~=0);
x = x(:,wh);
xx = xx(:,wh);

RES = [];

eb = pinv(x)*y;
er = eye(size(x,1))-(x*pinv(x));

df = [];
SS = [];
for ii = 1:length(F.err)
    if F.err(ii) > 0
        %%% Resize stuff if necessary
        a = F.FF{F.err(ii)}(Xind,:);
        a = a(:,find(nansum(a)>0));

        b =  [];
        ind5 = (ind2sub_aps(size(F.FF),F.err(ii)));
        ind5 = ind5([1 1+find(diff(ind5)>0)]);
        for jj = 1:length(ind5);
           tmp = F.FF{ind5(jj),ind5(jj)}(Xind,:);
           b(jj) = numel(find(nansum(tmp)>0));           
        end
        
        [r c cc] = MakeContrastMatrix('a',a,b,xx,x);
        
        eCons{ii} = c;
        
        mxm = (r-er);
        df(ii,1) = size(c,2)-numel(find(F.CovarCols)==1);
        
        
        for jj = 1:size(F.adj,2)
            if(F.adj(ii,jj) == 0)
                continue;
            end
            
            a = F.FF{F.adj(ii,jj)}(Xind,:);
            a = a(:,find(nansum(a)>0));
            
            b = F.adjLevs{ii,jj};

            [r c2 cc2] = MakeContrastMatrix('a',a,b,xx,x);
            mxm = mxm-(r-er);
            df(ii,1) = df(ii,1)-size(c2,2);
        end
        
        if nargout == 2 || F.writeRes == 1;
            [SS{ii,1} RES{ii}] = LoopEstimate(eb,x,mxm);
            RES{ii} = RES{ii}./repmat(sqrt(SS{ii}/df(ii))',1,size(RES{ii},2));            
            if F.writeRes==1
                fprintf('\n');
                
                tm = [pwd '/ResAll_00'];
                n = num2str(ii);
                tm((end-length(n))+1:end) = n;
                
                warning off;
                delete([tm '.nii']); delete([tm '.mat']);
                warning on;
                
                v = F.v;
                v.fname = [tm '.nii'];
                
                vol = nan(v.dim);
                for jj = 1:size(RES{ii},2);
                    vol(F.vec) = RES{ii}(:,jj);
                    v.n = [jj 1];
                    spm_write_vol(v,vol);
                end
                clear RES;

                h = spm_vol(v.fname);
                [FWHM,VRpv,R] = spm_est_smoothness(h,spm_vol('AllMask.nii'),[numel(h) df(ii)]);
                F.FWHM{ii} = FWHM;
                %F.ReslInfo{ii} = R;
                [tmpM tmpH] = openIMG('RPV.img');
                tmpH.fname = ['RPV' tm(end-2:end) '.nii'];
                spm_write_vol(tmpH,tmpM);
                delete('RPV.*')
            end
        else
            SS{ii,1} = LoopEstimate(eb,x,mxm);
        end
    elseif F.err(ii) == -1
        er = eye(size(x,1))-(x*pinv(x));
        
        if size(x,2) == 1
            df(end+1) = size(x,1)-1;
            eCons{ii} = 'Residual forming matrix';
        else
            tx = xx; tx(isnan(tx))=0;
            [z1 z2 z3] = svd(tx);
            tol = max(size(x))*max(abs(diag(z2)))*eps;
            df1 = sum(diag(z2)>tol);
            df1 = size(x,1)-df1;
            df(end+1) = df1;
            eCons{ii} = 'Residual forming matrix';
        end
        
        if nargout == 2 || F.writeRes == 1;
            [SS{ii,1}  RES{ii}] = LoopEstimate(y,1,er);
            RES{ii} = RES{ii}./repmat(sqrt(SS{ii}/df(ii))',1,size(RES{ii},2));
            %keyboard
            if F.writeRes==1
                %%%
                fprintf('\n');
                tm = [pwd '/ResAll_00'];
                n = num2str(ii);
                tm((end-length(n))+1:end) = n;
                warning off;
                delete([tm '.nii']); delete([tm '.mat']);
                warning on;
                v = F.v;
                v.fname = [tm '.nii'];
                vol = nan(v.dim);
                for jj = 1:size(RES{ii},2);
                    vol(F.vec) = RES{ii}(:,jj);
                    v.n = [jj 1];
                    spm_write_vol(v,vol);
                end
                %%%
                clear RES;

                h = spm_vol([tm '.nii']);
                [FWHM,VRpv,R] = spm_est_smoothness(h,spm_vol('AllMask.nii'),[numel(h) df(ii)]);
                F.FWHM{ii} = FWHM;
                %F.ReslInfo{ii} = R;
                [tmpM tmpH] = openIMG('RPV.img');
                tmpH.fname = ['RPV' tm(end-2:end) '.nii'];
                spm_write_vol(tmpH,tmpM);
                delete('RPV.*')
            end
        else
            SS{ii,1} = LoopEstimate(y,1,er);
        end
        
    elseif F.err(ii) == 0
        
    end
end

F.ErrorSS = SS;
F.ErrorDF = df;
F.eCons = eCons;

warning off
delete('AllMask.nii')
warning on