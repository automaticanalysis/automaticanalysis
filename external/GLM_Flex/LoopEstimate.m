function [out fullMat] = LoopEstimate(b,x,m)
%%% This computes the Sums of Squares.  Since there will typically be a
%%% lrage number of voxels it is most effecient to chunk the estimates
%%% rather than trying to estimate every voxel at once.  This script
%%% handles the chunking and SS computation.
%%%
%%% Inputs:
%%%
%%% b = the matrix of beta weights
%%% x = the design matrix (or whitened design matrix if that is being used)
%%% m = the mixing matrix.
%%%
%%% Outputs:
%%% out = a vector of Sums of Squares for contrast encoded in the mixing
%%% matrix m.
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

ints = 1:100:size(b,2);
sets = [];
sets(:,1) = ints(1:end);
sets(:,2) = [ints(2:end)-1 size(b,2)];

for ii = 1:size(sets,1)
    vec = sets(ii,1):sets(ii,2);
    out(vec) = diag((b(:,vec)'*x'*m*x*b(:,vec)));
end

if nargout==2
    fullMat = (b'*x'*m);
end