function cr = MatCorr(S1,S2);
%%% Get the pairwise correlations between two N by M data matrices.
%%%
%%% Inputs:
%%% S1 = a N observations x M Measurement data matrix
%%% S2 =  matrix of the same size as S1
%%%
%%% Outputs:
%%% cr = an M x 1 vector of the corelations between the paired columns M of
%%% S1 and S2.
%%%
%%% Example: 
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

m1 = repmat(mean(S1), size(S1,1),1);
m2 = repmat(mean(S2), size(S2,1),1);

S1 = S1-m1;
S2 = S2-m2;
out = LoopCross(S1,S2);
SS1 = LoopCross(S1,S1);
SS2 = LoopCross(S2,S2);
cr = (out./(sqrt(SS1.*SS2)))';