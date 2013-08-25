function [N newTarg] = SliceAndDice2(inputVol,Target,VoxDim,BB,interpOrd,prefix)
% function [N mat] = SliceAndDice(inputVol,Rot,VoxDim,BB,interpOrd,prefix)
%
% This function will also you to reslice a volume in a number of ways.  You
% can independently specify orientation, voxel size, and the bounding box
% for the output volume.  This script works by creating the desired affine
% (mat) matrix for the output volume.  Once this is in hand a direct
% mapping from input Image to out put image can be computed via
% Out.mat\In.mat, which is then passed to spm_slice_vol to do the actual
% reslicing.
%
% INPUTVOL:  The Volume to be rotated, resized, and resliced.  This can be
% a series so long as all images in the series have identical affine
% matrices.
% 
% ROT,VOXDIM,and BB can all be specied either numerically, with a target
% header (spm_vol), or a target filename.  Alternatively entering [], will
% use the INPUTVOL's own information to specify the parameters
% 
% ROT:  The rotation matrix or the target whose orientation you want to
% mimic. (4xF Euler/Affine matrix
%
% VOXDIM: The Voxel dimension to use for the image transformation.  These
% can be taken from the target, or specified manually (e.g. [3 3 3]).
%
% BB:  The Bounding Box to use for image transformations. Again this can be
% taken from the source, a target, or specified manually via a a 2x3 matrix
% (e.g. [-78 -112 -70; 78 76 90]
%
% INTERPORD: the interpolation order to use for resampling (see
% spm_slice_vol.m).  This is specified as an integer referencing the
% polynomial order to use for interpolation.  Additionally a second input
% specifying the fill value for padded voxels can also be specified.  For
% example 1 would use trilinear interpolation and [1 NaN] would use
% trilinear interploation, and any padding voxels added would be given a
% value of NaN.
%
% PREFIX: the filename profix to be applied when writing out the resliced
% image.  The default is 'SAD_' Slice and Dice!.
%
% Example:
% h1 = spm_vol('test1.img'); 
% h2 = spm_vol('test2.nii');
% [N mat] = SliceAndDice(h1,h2,h2,h2,[1 NaN]);
% spm_check_registration(spm_vol(char({h2.fname ['apsLice_' h1.fname]}))); shg
%
%
% Acknowledgements:  John Ashburner's reoirent.m script and various other
% scripts from John's Gems came in very handy in figuring out how to do
% this in a modular fashion.
%
% Written by Aaron P. Schultz June 3rd 2011
% aschultz@martinos.org
% Copyright (C) 2011,  Aaron P. Schultz
%
% Supported in part by the NIH funded Harvard Aging Brain Study (P01AG036694) and NIH R01-AG027435 
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

if nargin < 5;
    interpOrd = 0;
end

if nargin < 6;
    prefix = 'SAD_';
end

if isstruct(inputVol)
    h1 = inputVol;
else
    h1 = spm_vol(inputVol);
end
if numel(h1)>1
    tmp = [h1.mat]';
    if mean(mean(abs(diff([tmp(1:4:end,:) tmp(2:4:end,:) tmp(3:4:end,:) tmp(4:4:end,:)])))) ~= 0
       error('Volume Orientations are not all the same. Each Volume must be entered separately'); 
    end
end


if isempty(Target)
    TM = h1(1).mat;
    h2 = h1(1);
else
    if ischar(Target)
        h2 = spm_vol(Rot); h2 = h2(1);
    else
        h2 = Target;
    end
    TM = h2(1).mat;
end
tmp = spm_imatrix(TM);
TM1 = spm_matrix([0 0 0 tmp(4:6) sign(tmp(7:9))]);

if isempty(VoxDim)
    h3 = h1(1);
    tmp = spm_imatrix(TM1\h1(1).mat);
    voxdim = spm_imatrix(h1(1).mat*inv(spm_matrix([0 0 0 tmp(4:6) sign(tmp(7:9))])));
    voxdim = abs(voxdim(7:9));
else
    if isnumeric(VoxDim)
        voxdim = abs(VoxDim);
    end
    if ischar(VoxDim)
        h3 = spm_imatrix(VoxDim);
        tmp = spm_imatrix(TM1\h3(1).mat);
        voxdim = spm_imatrix(h3(1).mat*inv(spm_matrix([0 0 0 tmp(4:6) sign(tmp(7:9))])));
        voxdim = abs(voxdim(7:9));
    end
    if isstruct(VoxDim)
        h3 = VoxDim;
        tmp = spm_imatrix(TM1\h3(1).mat);
        voxdim = spm_imatrix(h3(1).mat*inv(spm_matrix([0 0 0 tmp(4:6) sign(tmp(7:9))])));
        voxdim = abs(voxdim(7:9));
    end
end

tmp = spm_imatrix(TM1);
tmp2 = spm_imatrix(TM);
TM2 = spm_matrix([0 0 0 tmp(4:6) voxdim.*sign(tmp2(7:9))]);
if isempty(BB);
    [bb bbb] = world_bb(h1(1));
    %%%  I need to find a way to pull the right values out of bbb
else
    if isnumeric(BB)
        bb = BB;
%         [0 0 0 1; 0 0 1 1; 0 1 0 1; 0 1 1 1; 1 0 0 1; 1 0 1 1; 1 1 0 1; 1 1 1 1]*TM1'
    end
    if ischar(BB)
        h4 = spm_vol(BB);
        [bb bb2] = world_bb(h4(1));
    end
    if isstruct(BB)
        [bb bb2] = world_bb(BB);
    end
end
tb = world_bb(h2);
tmp = h2.mat;
tmp(1:3,4) = 0;
tar = h2.mat(:,4)'+([1 1 1 1]*tmp');

% try
origin = [];
for ii = 1:3;
    [r c] = find(tb==tar(ii));
    origin(ii) = bb(r(1),c(1));
end
% catch
%     keyboard; 
% end
origin = [origin 1] - ([1 1 1 1]*TM2');
origin = origin(1:3);
    
tmp = spm_imatrix(TM2);
TM3 = spm_matrix([origin tmp(4:6) voxdim.*sign(tmp2(7:9))]);
dims = round(abs(diff([bb [1;1]]*inv(TM3'))))+1;

newTarg = TM3;
imgdim = dims;
tm = newTarg\h1(1).mat;

%%% Use spm_slice_vol to reslice the image to the desired specifications
NN = zeros([imgdim numel(h1)]);
for jj = 1:length(h1);
    
    V = h1(jj);
    V = rmfield(V,'private');
    V.fname = strtrim([prefix V.fname]);
    V.dim = imgdim(1:3);
    V.mat = newTarg;

    N = zeros(imgdim);
    for i = 1:imgdim(3) 
        M = tm\spm_matrix([0 0 i]);
        img = spm_slice_vol(h1(jj), M, imgdim(1:2), interpOrd);
        N(:,:,i) = img;
    end
    NN(:,:,:,jj) = N;
    if ~isempty(prefix)
        spm_write_vol(V,N);
    end
end

if numel(NN) > numel(N)
    N = NN;
end

return
%%
clear
cd /autofs/space/plato_004/users/RestingStateScans/R01_Long/RLB008C4_MB_3m;
struc = getNIFTI('RLB008C4_MB_3m','/MPRAGE_1.nii','task')

h1 = spm_vol(struc{1});
h2 = spm_vol('meanst_RLB008C4_MB_3m_Rest.nii');
% M  = spm_matrix(x1); 

% return
[N mat] = SliceAndDice3(h2,h1,h2,[],1,'ac');
close all; OrthoView({{struc{1}} {['ac' 'meanst_RLB008C4_MB_3m_Rest.nii']}});

