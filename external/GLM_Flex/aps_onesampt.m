function aps_onesampt(nn, minN, RO, ROT)
%%% Perform a one-sample t-test using GLM_Flex.
%%%
%%% Inputs
%%% nn = a cell array of filenames.  
%%% 
%%% Optional Inputs:
%%% minN =  the minimum number of observations per condition needed to 
%%% analyze a voxel.  Default = All data must be present.
%%%
%%% RO = Flag to remove outliers 1 = yes, 2 = no.  Default = no
%%%
%%% ROT = outlier threshold.  Default this will be set automatically based
%%% on N size and design if RO is turned on.  Setting this yourself is not
%%% reccomended.
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

if nargin < 2
    minN = numel(nn);
    RO = 0;
elseif nargin < 3
    RO = 0;
end


IN.N_subs = [numel(nn)];
IN.Between = [1];
IN.EqualVar = [1];
IN.Independent = [1];

IN.Scans = nn;
F = CreateDesign(IN);

I.OutputDir = pwd;
I.F = F;
I.Scans = nn;
I.minN = minN;
I.RemoveOutliers = RO;
if nargin == 4
    I.Thresh = ROT;
end

I.Cons.Groups = {1};
I.Cons.Levs = 1;
I.Cons.mean = 1;
I.Cons.ET = 1;

I = GLM_Flex(I);


% I = GLM_Flex_Contrasts(I);
