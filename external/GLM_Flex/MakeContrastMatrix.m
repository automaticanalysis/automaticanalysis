function [r c cc rr] = MakeContrastMatrix(type,sDM,levs,xx,x,CovarCols)
%%% Generate contrast matricies and mixing matrices
%%%
%%% Inputs:
%%%
%%% type = the method of generating the pre contrasts. In general it is
%%% advisable to always use 'a', this seems to be the most general and
%%% stable method.  The other methods will likely be removed in future
%%% versions
%%%
%%% sDM = Sub Design Matrix. here you pass in the columns of the design
%%% matrix corresponding to the effect that you are interested in (each
%%% condition should be uniquely specified by a column).
%%%
%%% levs = the number and distribution of levels.  A four level on way
%%% anova would be specified as [4] wheras a 2x2 interation would be
%%% specified as [2 2].  The order should be the same as in the model. That
%%% is if you put in the 3 level factor before a 2 level factor then the
%%% levs order would be [3 2] and not [2 3].  Order here implies the order
%%% of looping to create the conditions "for ii = 1:3; for jj = 1:2; ..."
%%% will produce a different design matrix patternd than "for ii = 1:2; for
%%% jj = 1:3; ..." so be sure that the order is correct or the contrast
%%% will not be correct.
%%%
%%% xx = the full design matrix with NaN's in place of zeros.
%%%
%%% x = the whitened designed matrix or if no variance correction, it will
%%% be the same as xx only with zeros instead of NaN's
%%%
%%% CovarCols = a vector of ones and zeros coding whether or not each
%%% column of the design matrix is a covariate.  If not specified it is
%%% assumed that there are no covaraites or that covariates should be
%%% included in the group/condition definitions (e.g. don't remove the
%%% effect of covariate from group).
%%%
%%% Outputs:
%%%
%%% r = the mixing matrix to compute the SS for the contrast. This is used
%%% in SS = b'*x'*m*x*b.  m is the mixing matrix which is r-rr.
%%%
%%% rr =  the residual forming matrix for the full model
%%%
%%% c = the contrasts vector/matrix
%%%
%%% cc = the pre-contrast vector/matrix.  
%%%
%%% Note: c is comuted from cc.  Based on the number of levels and their
%%% distribution across factors, the contrast c is created by folding cc
%%% and the differentiating across the proper dimensinos
%%% 
%%%
%%% Written by Aaron Schultz, May 5th, 2010 (aschultz@martinos.org)
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


if numel(isnan(xx))==0;
    xx(find(xx==0))=NaN;
end

if nargin<5 || isempty(x)
    x = xx;
    x(isnan(x)) = 0;
end

if nargin<6 || isempty(CovarCols)
    CovarCols = zeros(1,size(x,2));
end

levs = levs(end:-1:1);

if     type == 'a';  %% Exact matches-covariates are not included in factor contrasts
    cc = [];
    for ii = 1:size(sDM,2)
        ind1 = find(~isnan(sDM(:,ii)));
        go = 0;
        for jj = 1:size(xx,2);
            ind2 = find(~isnan(xx(:,jj)));
            ind = intersect(ind1,ind2);
            if isempty(ind)
                go = 0;
            else
                match = mean(sDM(ind,ii)==xx(ind,jj));
                if match == 1;
                    go = 1;
                else
                    go = 0;
                end
            end
            if go == 1
                if CovarCols(jj)==0
                    cc(jj,ii) = sum(xx(ind,jj))./numel(ind1);
                else
                    cc(jj,ii) = sum(~isnan(xx(ind,jj)))./numel(ind1);
                end
            else
                cc(jj,ii) = 0;
            end
        end
    end
elseif type == 'b';  %% Get mean of overlapping columns
    c = [];
    for ii = 1:size(sDM,2)
        ind = find(~isnan(sDM(:,ii)));
        tmp = (xx(ind,:));
        tmp(isnan(tmp))=0;
        cc(:,ii) = mean(tmp);
    end
elseif type == 'c';
    c = [];
    for ii = 1:size(sDM,2)
        ind = find(~isnan(sDM(:,ii)));
        tmp = (xx(ind,:));
        tmp(:,find(CovarCols==1)) = ~isnan(tmp(:,find(CovarCols==1)));
        tmp(isnan(tmp))=0;
        cc(:,ii) = mean(tmp);
    end
elseif type == 'd';
    
elseif type == 'e';
    
end


tmp = reshape(cc,[size(cc,1) levs]);
tmp = squeeze(tmp);
if size(tmp,2)==1
    c = tmp;
else
    c = differencer(tmp)*-1;
%     if size(c,2)==1;
%         c = c*-1;
%     end
end
    
c0 = eye(size(x,2))-(c*pinv(c));
x0 = x*c0;
r = eye(size(x0,1))-(x0*pinv(x0));

c0 = eye(size(x,2))-(cc*pinv(cc));
x0 = x*c0;
rr = eye(size(x0,1))-(x0*pinv(x0));

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
