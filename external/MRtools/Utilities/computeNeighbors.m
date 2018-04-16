function A = computeNeighbors(im1,im2)
%%% Compute the number of Neihbors for each voxel in the mask im1 and then
%%% return the number of neighbors for in im1 for each voxel in mask im2
%%% (note that this involves conjoing the im1 and im2)
%%%
%%% Example A = computeNeighbors('/autofs/space/madrc_007/users/FreeSurfer/Mixed_Design_young_struc/Yng_mixed_10/label/wbil.cortexhip_native_maskoutsus.img','/space/madrc/7/users/FreeSurfer/R01_structural/RLB001C4_MA_3m/label/pcc5mmsph_0_-53_26.img')
%%%
% im1 = '/autofs/space/madrc_007/users/FreeSurfer/Mixed_Design_young_struc/Yng_mixed_10/label/wbil.cortexhip_native_maskoutsus.img';
% im2 = '/space/madrc/7/users/FreeSurfer/R01_structural/RLB001C4_MA_3m/label/pcc5mmsph_0_-53_26.img';
% A = computeNeighbors(im1,im2)
%%% Written by Aaron Schultz (aschultz@martinos.org)
%%% Copyright (C) 2014,  Aaron P. Schultz
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

m =  openIMG(im1);
m2 = openIMG(im2);

V = square2vect3D(2:size(m,1)-1,2:size(m,2)-1,2:size(m,3)-1,size(m));
% ind = find(m(V)~=0);

%%% In Plane 
comp = [m(V-1) m(V+1) m(V-size(m,1)) m(V-size(m,1)+1) m(V-size(m,1)-1) m(V+size(m,1)) m(V+size(m,1)+1) m(V+size(m,1)-1)];
%%% In back of plane.
B = (V-(size(m,1)*size(m,2)));
comp = [comp m(B) m(B-1) m(B+1) m(B-size(m,1)) m(B-size(m,1)+1) m(B-size(m,1)-1) m(B+size(m,1)) m(B+size(m,1)+1) m(B+size(m,1)-1)];
%%% In front of plane.
F = (V+(size(m,1)*size(m,2)));
comp = [comp m(F) m(F-1) m(F+1) m(F-size(m,1)) m(F-size(m,1)+1) m(F-size(m,1)-1) m(F+size(m,1)) m(F+size(m,1)+1) m(F+size(m,1)-1)];

% K1 = sum(comp(ind,:),2)+1;
T = zeros(size(m));
T(V) = sum(comp,2)+1;



A.VoxInRoi = numel(find(m2>0));
A.VoxInConj = numel(find((m.*m2)>0));
A.K_values = T(find((m.*m2)>0));
A.Index1 = find((m.*m2)>0);
A.Index2 = vect2square3D(A.Index1,size(m));

% aa = find((m1.*m2)>0);
% aa = find(m2>0);
% ch = openIMG('/autofs/space/plato_001/users/liang/Young_mixed_reho/Yng_mixed_10/structurals/resting_state/SPM_stats_filter_noSM/NUIROIreg_bwWM_3/test_K_corticalhipmask/reho_Ksize.img');
% [ch(aa) T(aa)]

