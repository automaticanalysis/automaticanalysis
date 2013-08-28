function aps_pairedt(nn1, nn2, minN,RO,ROT)
%%% Perform a paired-samples t-test using GLM_Flex
%%%
%%% Inputs
%%% nn1 = cell array of filenames for sample #1
%%% nn2 = cell array of filenames for sample #2
%%%     Note: the oder of subjects in nn1 and nn2 must be the same
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

if nargin < 3
    minN = numel(nn1);
    RO = 0;
    ROT = [];
elseif nargin < 4
    RO = 0;
    ROT = [];
elseif nargin < 5
    ROT = [];
end
IN.N_subs = [numel(nn1)];
IN.Within = [2];
IN.EqualVar = [1];
IN.Independent = [1];

nn = interleave(nn1,nn2);

F = CreateDesign(IN);

I.OutputDir = pwd;
I.F = F;
I.Scans = [nn];
I.minN = minN;
I.RemoveOutliers = RO;
I.Thresh = ROT;
I.CompOpt=0;

% I.writeRes=0;

I.Cons(1).Groups = {1 2};
I.Cons(1).Levs = 2;
I.Cons(1).mean = 1;
I.Cons(1).ET = 2;

I.Cons(2).Groups = {1};
I.Cons(2).Levs = 1;
I.Cons(2).mean = 0;
I.Cons(2).ET = 1;

I.Cons(3).Groups = {2};
I.Cons(3).Levs = 1;
I.Cons(3).mean = 0;
I.Cons(3).ET = 1;

I = GLM_Flex(I);

%I = GLM_Flex_Contrasts(I);


