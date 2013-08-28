function [vals o1 o2 flips] = matInfo(mat)
%%% Won't go into too much detail on this one.  This get's some information
%%% from the affine matrix in the image headers.  I use this info in the
%%% SliceAndDice Script.
%%%
%%% Input:
%%% mat = affine matrix from volume header
%%%     e.g. h = spm_vol('Image.img'); mat = h.mat;
%%%
%%% Outputs:
%%% Won't go into detail on these at the moment.  This script is primarily
%%% used by SliceAndDice.m
%%%
%%% Written by Aaron Schultz - aschultz@martinos.org
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
% keyboard;
tmp = spm_imatrix(mat);
tmp = tmp(7:9);
tmp = round(tmp*10000)/10000;
oris = [[mat(1,1) mat(2,2) mat(3,3)];...
        [mat(1,1) mat(3,2) mat(2,3)];...
        [mat(2,1) mat(1,2) mat(3,3)];...
        [mat(2,1) mat(3,2) mat(1,3)];...
        [mat(3,1) mat(1,2) mat(2,3)];...
        [mat(3,1) mat(2,2) mat(1,3)]];
ords1 = [1 2 3;
         1 3 2;
         2 1 3;
         3 1 2;
         2 3 1;
         3 2 1];
ords2 = [1 2 3;
         1 3 2;
         2 1 3;
         2 3 1;
         3 1 2;
         3 2 1];


oris = round(oris*10000)/10000;     


match = (oris)==(repmat(tmp,6,1));
ori = find(mean(match,2)==1);

if isempty(ori)
    match = abs((oris))==abs(repmat((tmp),6,1));
    ori = find(mean(match,2)==1);
end

if isempty(ori)
    ori = find(sum(abs(abs(oris)-abs(repmat(tmp,6,1))),2)==min(sum(abs(abs(oris)-abs(repmat(tmp,6,1))),2)));
end


vals = oris(ori,:);
o1 = ords1(ori,:);
o2 = ords2(ori,:);

% flips = sign(oris(ori,:))-[-1 1  1];
flips = sign(oris(ori,:));