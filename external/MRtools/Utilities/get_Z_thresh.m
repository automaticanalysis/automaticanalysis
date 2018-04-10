function [zthresh r2 r] = get_Z_thresh(N,NR,SigLevel)
% N = 20; %Sample Size (Number of Subjects for second level analysis, number of time points for first level analysis
% NR = 2; %Number of Regressors 2 for a Seed based analysis the seed plus a constant in the regressions equation.
% SigLevel = .001;
% N = length(SPM.xY.P);
% NR = size(SPM.xX.X,2);
%%% Written by Aaron Schultz (aschultz@martinos.org)
%%%
%%% Copyright (C) 2012,  Aaron P. Schultz
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
 
df1 = NR-1;
df2 = N-df1-1;

Fs = 0:.001:100;
F_val = 1-cdf('F',Fs,df1,df2);
ind = find(F_val<SigLevel);
minF = Fs(ind(1));

 
r2 = 1/((1/(minF/df2*df1))+1);
r = sqrt(r2);
zthresh = atanh(r);
% [r2 r zthresh]

 