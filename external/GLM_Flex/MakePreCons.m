function cc = MakePreCons(sDM,xx,CovarCols)
%%% This os a fucntion to get vectors and matricies for computing means
%%% from beta weights in GLMs
%%%
%%% Inputs:
%%% sDM = Sub Design Matrix. here you pass in the columns of the design
%%% matrix corresponding to the effect that you are interested in (each
%%% condition should be uniquely specified by a column).
%%%
%%% xx = the full design matrix with NaN's in place of zeros.
%%%
%%% CovarCols = a vector of ones and zeros coding whether or not each
%%% column of the design matrix is a covariate.  If not specified it is
%%% assumed that there are no covaraites or that covariates should be
%%% included in the group/condition definitions (e.g. don't remove the
%%% effect of covariate from group).
%%%
%%% Outputs:
%%% cc = a matrix with each column containing the weights for each column
%%% of the design matrix for that particular condition.  If sDM has four
%%% columns cc will have four columns.  The best way to think of cc is as a
%%% way of recovering group means from the betas.  cc'*betas will give you
%%% the group means.
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

if numel(find(isnan(xx)))==0;
    xx(find(xx==0))=NaN;
end

if numel(find(isnan(sDM)))==0;
    sDM(find(sDM==0))=NaN;
end

if nargin<3 || isempty(CovarCols)
    CovarCols = zeros(1,size(xx,2));
end

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
