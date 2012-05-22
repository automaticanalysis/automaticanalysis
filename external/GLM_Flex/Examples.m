%%% Examples on how to Create GLM designs with CreateDesign.m
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

%% 2x2 Between ANOVA
IN.N_subs = [5 5 5 5];
IN.Between = [2 2];
IN.Interactions = {[1 2]};
IN.EqualVar = [0 0];
IN.Independent = [1 1];

F = CreateDesign(IN);

figure(100); clf;
X = F.XX;
X(isnan(X))=0;
imagesc(X); colormap(gray); shg

%% 2x2 Within ANOVA
clear IN F
IN.N_subs = [10];
IN.Within = [2 2];
IN.Interactions = {[1 2]};
IN.EqualVar = [0 0];
IN.Independent = [0 0];

F = CreateDesign(IN);

figure(100); clf;
X = F.XX;
X(isnan(X))=0;
imagesc(X); colormap(gray); shg

%% 2x3x2 ANOVA
clear IN F
IN.N_subs = [5 5 5 5 5 5 5 5 5 5 5 5];
IN.Between = [2 3 2];
IN.Interactions = {[1 2] [1 3] [2 3] [1 2 3]};
IN.EqualVar = [0 0 0];
IN.Independent = [1 1 1];

F = CreateDesign(IN);

figure(100); clf;
X = F.XX;
X(isnan(X))=0;
imagesc(X); colormap(gray); shg

%% 2x2 Mixed ANOVA
clear IN F;
IN.N_subs = [5 5];
IN.Between = [2];
% IN.BetweenLabs = {{'A' 'B'}};
IN.Within = [2];
% IN.WithinLabs = {{'C' 'D'}};
% IN.FactorLabs = {'Groups' 'Pre/Post'};
IN.Interactions = {[1 2]};
IN.EqualVar = [0 0];
IN.Independent = [1 1];

F = CreateDesign(IN);

figure(100); clf;
X = F.XX;
X(isnan(X))=0;
imagesc(X); colormap(gray); shg
%% One between and two Within Factors
clear IN F;
IN.N_subs = [5 5];
IN.Between = [2];
IN.Within = [2 2];
IN.Interactions = {[1 2] [1 3] [2 3] [1 2 3]};
IN.EqualVar = [1 1 1];
IN.Independent = [1 1 1];

F = CreateDesign(IN);

figure(100); clf;
X = F.XX;
X(isnan(X))=0;
imagesc(X); colormap(gray); shg
%%
xx = F.XX;
% [r c cc rr] = MakeContrastMatrix('a',F.FF{3,3},4,xx,xcc2);
cc2 = MakePreCons([[ones(1,11) zeros(1,9)]' [zeros(1,11) ones(1,9)]'],F.XX)
cc3 = MakePreCons(F.FF{1,1},F.XX)
