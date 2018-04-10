function N = resizeVol(V1,V2)
%%% This function is based on the spm_mask.m function.  V1 is the volume
%%% being resized, V2 is the target Volume being resized to;
%%% 
%%% Note: The inputs are the structures returned from spm_vols.
%%%
%%% Written by Aaron Schultz (aschultz@martinos.org)
%%%
%%% Last Updated Dec. 11 2012;
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

M   = V2(1).mat;
dim = V2(1).dim(1:3);

N = zeros(V2(1).dim);
for j=1:dim(3)
    msk = zeros(dim(1:2));
    Mi  = spm_matrix([0 0 j]);
    M1  = M\V1.mat\Mi;
    img = spm_slice_vol(V1,M1,dim(1:2),[0 NaN]);    
    msk = msk + (img~=0 & isfinite(img));
%     msk = msk+img;
    N(:,:,j) = msk;
end

N(N==0) = NaN;

