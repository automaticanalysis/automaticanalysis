function aps_onewayanova_repeated(nn,minN,RO,ROT)
%%% Perform a one-way repeated measures ANOVA using GLM_Flex
%%%
%%% Inputs
%%% nn = double cell array of filenames.  The number of cells in the outer
%%% array will be equal to the number of levels.  Each inner array should
%%% contain a list of file names for that condition.  The order of subjects
%%% needs to be the same in each inner array:  Example:  
%%% { {'S1C1.img' 'S2C1.img' ...} {'S1C2.img' 'S2C2.img' ...} }
%%%     Note: the oder of subjects in each cell must be the same
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
    minN = [];
    for ii = 1:length(nn)
        minN = [minN numel(nn{ii})];
    end
    minN = min(minN);
    RO = 0;
    ROT = 0;
elseif nargin < 3
    RO = 0;
    ROT = 0;
elseif nargin < 4
    ROT = [];
end

Scans = [];
IN.N_subs  = numel(nn{1});
com = 'Scans = interleave(';
for ii = 1:length(nn)
    com = [com 'nn{' num2str(ii) '},'];
end
com(end:end+1) = ');';
eval(com);

IN.Within = [numel(nn)];
IN.EqualVar = [1];
IN.Independent = [1];

F = CreateDesign(IN);

I.OutputDir = pwd;
I.F = F;
I.Scans = Scans;
I.minN = minN;
I.RemoveOutliers = RO;
I.Thresh = ROT;

I = GLM_Flex(I);

g = [];
for ii = 1:length(nn)
    g{ii} = ii;
    
    I.Cons(ii+1).Groups = {ii};
    I.Cons(ii+1).Levs = 1;
    I.Cons(ii+1).mean = 0;
    I.Cons(ii+1).ET = 1;
end

I.Cons(1).Groups = g;
I.Cons(1).Levs = numel(g);
I.Cons(1).mean = 1;
I.Cons(1).ET = 2;

I = GLM_Flex_Contrasts(I);

