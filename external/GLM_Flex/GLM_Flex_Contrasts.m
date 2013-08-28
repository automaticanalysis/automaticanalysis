function [I oo] = GLM_Flex_contrasts(I,these,oo)
%%% This script generates the contrast images for the GLM_Flex package.
%%% Go to: http://nmr.mgh.harvard.edu/harvardagingbrain/People/AaronSchultz/Aarons_Scripts.html
%%% for more information on this script and how to use it.
%%%
%%% After running I = GLM_Flex(I), add the Cons structure to I and
%%% then pass it into this scrtip.  The Cons field can be specified as one
%%% of the following:
%%%
%%% I.Cons(i).name = '';
%%% I.Cons(i).Groups = [];
%%% I.Cons(i).Levs = [];
%%% I.Cons(i).ET = [];
%%% I.Cons(i).mean = 0;
%%%
%%% OR
%%%
%%% I.Cons(i).name = '';
%%% I.Cons(i).c = [];
%%% I.Cons(i).ET = [];
%%% I.Cons(i).mean = [];
%%%
%%% Copyright (C) 2011,  Aaron P. Schultz
%%% Written by Aaron Schultz (aschultz@martinos.org)
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


Cons = I.Cons;

if nargin==1
    these = 1:numel(Cons);
end

MeanVol = [];

cd(I.OutputDir);
if nargin<3
    load([pwd filesep 'oo.mat']);
end
for ii = these;
    lab = '0000';
    f = num2str(ii);
    lab(end-length(f)+1:end)=f;
    delete(['*' lab '*.nii'])
    
    if ~isfield(Cons(ii),'mean');
        Cons(ii).mean=0;
    end
    if ~isfield(Cons(ii),'name') || isempty(Cons(ii).name);
        Cons(ii).name=['Con'];
    end
    if ~isfield(Cons(ii),'ET');
        Cons(ii).ET=[];
    end
    
    tx = I.F.XX;
    txx = I.aX;
    
    if ~isfield(Cons(ii),'c') || isempty(Cons(ii).c)
        c = getC(Cons(ii),tx,I.F.CovarCols);
    else
        c = Cons(ii).c;
    end
    c0 = eye(size(txx,2))-(c*pinv(c));
    x0 = txx*c0;
    r = eye(size(x0,1))-(x0*pinv(x0));
    
    fprintf('\n\n');
    if ~isempty(Cons(ii).ET);
        errLev = Cons(ii).ET;
        disp(['Using ResMS #' num2str(errLev)]);
        df1 = size(c,2);
        df2 = oo.EDF{1}(errLev);
    else
        if I.F.err(1) == -1;
            errLev = 1;
            disp(['Using ResMS #' num2str(errLev)]);
            df1 = size(c,2);
            df2 = oo.EDF{1}(errLev);
        else
            errLev = GetErrLev(c,I.F);
            if isempty(errLev);
                errLev = input('Enter which error term to use:   ');
            end
            disp(['Using ResMS #' num2str(errLev)]);
            df1 = size(c,2);
            df2 = oo.EDF{1}(errLev);
        end
    end
    I.Cons(ii).ET = errLev;
    h = I.v;
    Vol1 = zeros(h.dim)*NaN;
    Vol2 = zeros(h.dim)*NaN;
    Vol3 = zeros(h.dim)*NaN;
    Vol4 = zeros(h.dim)*NaN;
    
    
    h.n = [1 1];
    h.dt = [16 0];
    if size(c,2)>1 || numel(Cons(ii).Levs)>1
        h.descrip = ['SPM{F_' '[' num2str(size(c,2)) ',' num2str(oo.EDF{1}(errLev)) ']} - computed with APS Estimate'];
    else
        h.descrip = ['SPM{T_' '[' num2str(df2) ']} - computed with APS Estimate'];
    end
    
    if Cons(ii).mean==1;
        for kk = 1:prod(Cons(ii).Levs);
            MeanVol{ii}{kk} = Vol3;
        end
    end
    
    persisText
    for jj = 1:length(oo.Xind)
        if mod(jj,100)==0
            persisText(['Processing Model ' num2str(jj) ' of ' num2str(length(oo.Xind))]);
        end
        
        Xind = oo.Xind{jj};
        beta = oo.beta1{jj};
        
        x = I.F.XX(Xind,:);
        xx = oo.W{jj}*oo.X{jj};
        
        if ~isfield(Cons(ii),'c') || isempty(Cons(ii).c)
            [c cc] = getC(Cons(ii),x,I.F.CovarCols);
        else
            c = Cons(ii).c;
            cc = [];
        end
        
        c0 = eye(size(xx,2))-(c*pinv(c));
        x0 = xx*c0;
        r = eye(size(x0,1))-(x0*pinv(x0));
        
        rr = eye(size(xx,1))-(xx*pinv(xx));
        
        try
            oo.Cons(ii).c{jj} = c;
        catch
            oo = rmfield(oo,'Cons');
            oo.Cons(ii).c{jj} = c;
        end
        
        SS = LoopEstimate(beta,xx,r-rr);

        Vol1(oo.vec{jj}) = SS./df1;
        try
            den = oo.ESS{jj}{errLev}./oo.EDF{jj}(errLev);
        catch; continue; end;
        
        F1 = ((SS./df1)./den);
        
        df11 = df1;
        df21 = oo.EDF{jj}(errLev);
        
        if ~isreal(F1)
            for zz = 1:numel(F1)
               if ~isreal(F1(zz))
                   F1(zz)=0;
               end
            end
        end
        
        try
            p1 = spm_Fcdf(F1,df11,df21);
            F = spm_invFcdf(p1,df1,df2);
            F(find(~isfinite(F))) = F1(find(~isfinite(F)));
            
            if size(c,2)>1 || numel(Cons(ii).Levs)>1
                Vol2(oo.vec{jj}) = F;
                filename = [lab '_F_' Cons(ii).name '.nii'];
            else
                Vol2(oo.vec{jj}) = sqrt(F).*sign(c'*beta);
                Vol4(oo.vec{jj}) = c'*beta;
                filename = [lab '_T_' Cons(ii).name '.nii'];
            end
            
        catch
            disp('err');
        end
        
        if Cons(ii).mean == 1;
            if exist('cc') ~= 0
                for kk = 1:size(cc,2);
                    MeanVol{ii}{kk}(oo.vec{jj}) = cc(:,kk)'*beta;
                end
            end
        end
    end
 
    fprintf('\n\n');
    writeIMG(h,Vol2,filename);

    
    
    filename = ['ess_' lab '.nii'];
    h.descrip = [];
    writeIMG(h,Vol1,filename);
    
    for jj = 1:length(MeanVol)
        for kk = 1:length(MeanVol{jj});
            filename = ['MeanVol_C' num2str(jj) '_G' num2str(kk) '.nii'];
            writeIMG(h,MeanVol{jj}{kk},filename);
        end
    end
    
    if size(c,2) == 1;
        filename = ['con_' lab '.nii'];
        writeIMG(h,Vol4,filename);
    end
end

if I.writeI==1
    save I.mat I;
end

if nargin < 3
    save oo.mat oo;
end

end

function [c cc] = getC(Con,x,CovarCols)
    if iscell(Con.Groups)
        if iscell(Con.Groups{1})
            for zz = 1:length(Con.Groups)
                tmpX = zeros(size(x,1),1)*NaN;
                tmpX(Con.Groups{zz}{1})=1;
                cc = MakePreCons(tmpX,x,CovarCols);
                tc(:,zz) = mean(cc,2);
            end
        else
            tc = [];
            for zz = 1:length(Con.Groups)
                if size(Con.Groups{zz},1)>1
                    tmpX = Con.Groups{zz};
                else
                    tmpX = x(:,Con.Groups{zz});
                end
                cc = MakePreCons(tmpX,x,CovarCols);
                tc(:,zz) = mean(cc,2);
            end
        end

        levs = Con.Levs(end:-1:1);
        if levs == 0
            c = tc;
        else
            tmp = reshape(tc,[size(tc,1) levs]);
            tmp = squeeze(tmp);
            if size(tmp,2)==1
                c = tmp;
            else
                c = differencer(tmp);
                if numel(levs)==1
                    c = c*-1;
                end
            end
        end
        cc = tc;
    else
        [r,c,cc] = MakeContrastMatrix('a',Con.Groups,Con.Levs,x,[],CovarCols);
    end
end

function errLev = GetErrLev(c,F)

[row,col] = find(abs(c)==1);
ind = unique(row);

wh = F.name(ind);

tr = [];
for ii = 1:length(wh);
    tmp = wh{ii};
    in1 = find(tmp=='x');
    if isempty(in1);
        in2 = find(tmp=='_');
        fac = str2num(tmp(in2(1)+1:in2(2)-1));
    else
        in1 = [0 in1 numel(tmp)+1];
        fac = [];
        for jj = 1:length(in1)-1
            tmp2 = tmp(1+in1(jj):in1(jj+1)-1);
            in2 = find(tmp=='_');
            fac(jj) = str2num(tmp2(in2(1)+1:in2(2)-1));
        end
    end
    tr = [tr fac];
end
tr = unique(tr);

if numel(tr)==1;
    tr = [tr tr];
end

try
    ij = sub2ind_aps(size(F.FF),tr);
    [row,col] = find(F.eff == ij);
    errLev = row;
catch
    error('Cannot figure out which error term to use!');
end
end

function out = differencer(tmp,count)

if nargin == 1;
    count = numel(size(tmp));
end

out = diff(tmp,1,count);

if count ~= 2;
    out = differencer(out,count-1);
else
    if numel(size(out))>2
        ss = size(out);
        out = reshape(out,size(out,1),prod(ss(2:end)));
    end
    return
end
end