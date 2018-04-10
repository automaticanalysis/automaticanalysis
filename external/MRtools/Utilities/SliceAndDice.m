function [N mat] = SliceAndDice(inputVol,Rot,VoxDim,BB,interpOrd,prefix)
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

if isempty(Rot)
    rm = spm_matrix([0 0 0]);
else
    if isnumeric(Rot)
        rm = Rot;
    else
        if ischar(Rot)
            h2 = spm_vol(Rot); h2 = h2(1);
        else
            h2 = Rot;
        end
        
        tmp1 = spm_imatrix(h1(1).mat);
        tmp2 = spm_imatrix(h2(1).mat);
        [vals o1 o2 flips] = matInfo(h2.mat);
        i1 = sort(sub2ind([4 4], o1,o2));
        v1 = h2.mat(i1);
        h2.mat(i1) = sign(h2.mat(i1)).*abs(tmp1(7:9));
                
        mat = h2.mat\h1(1).mat;
        matp = h1(1).mat\h2.mat;
        
        [vals o1 o2 flips] = matInfo(mat);
        
%         mat(1:3,1:3) = sign(round(mat(1:3,1:3)./max(max(abs(mat(1:3,1:3))))));
        mat(1:3,4) = matp(o1(1:3),4);
        
        rm = mat;
    end
end

if isempty(VoxDim)
    %Made changes to the VoxDim Section on Sept. 20th 2011
%     tmp = spm_imatrix(h1(1).mat);
%     voxdim = abs(tmp(7:9));
%     [vals o1 o2 flips] = matInfo(h1.mat);
%     ord = o1;
    
      voxdim = spm_imatrix(h1.mat);
      voxdim = abs(voxdim(7:9));
      ord = 1:3;
      
%       h3 = h1(1);
%       tmp = spm_imatrix(TM\h1(1).mat);
%       voxdim = spm_imatrix(h1(1).mat*inv(spm_matrix([0 0 0 tmp(4:6)])));
%       voxdim = abs(voxdim(7:9));
else
    if isnumeric(VoxDim)
        voxdim = VoxDim;
        ord = 1:3;
    end
    if ischar(VoxDim)
%         h3 = spm_vol(VoxDim);
%         tmp = spm_imatrix(h3(1).mat);
%         voxdim = abs(tmp(7:9));
%         [vals o1 o2 flips] = matInfo(h3.mat);
%         ord = o1;

        voxdim = spm_imatrix(h3.mat);
        voxdim = abs(voxdim(7:9));
        ord = 1:3;
    end
    if isstruct(VoxDim)
%         %keyboard;  
%         tmp = spm_imatrix(VoxDim.mat);
%         voxdim = abs(tmp(7:9));
%         [vals o1 o2 flips] = matInfo(VoxDim.mat);
%         ord = o1;
        
        voxdim = spm_imatrix(VoxDim.mat);
        voxdim = abs(voxdim(7:9));
        ord = 1:3;
        
        [a b c d] = matInfo(h2.mat);
    end
end

if isempty(BB);
    [bb bb2] = world_bb(h1(1));
else
    if isnumeric(BB)
        bb = BB;
    end
    if ischar(BB)
        h4 = spm_vol(BB);
        [bb bb2] = world_bb(h4(1));
    end
    if isstruct(BB)
        [bb bb2] = world_bb(BB);
    end
end
%%% The below works by creating the output affine matrix.

%%% Affine Matrix after applying rotations
mat = h1(1).mat*inv(rm);

%%% Affine Matrix with desired voxel dimensions.
[vals o1 o2 flips] = matInfo(mat);
mat(sort(sub2ind([4 4],o1,o2))) = sign(mat(sort(sub2ind([4 4],o1,o2)))).*voxdim(ord);
% mat(sort(sub2ind([4 4],o1,o2))) = sign(mat(sort(sub2ind([4 4],o1,o2)))).*voxdim;

%%% Affine Matrix with desired bounding box
tmp = sign([1 1 1 1]*mat'); 
ib = find(sum(abs(sign(bb2')-repmat(tmp(1:3),8,1)),2)==0);
tmp = bb2(:,ib)';
tmp = tmp+(sign(tmp).*voxdim(ord));
mat(1:3,4) = tmp;

% keyboard;
%%% Figure out the output image matrix dimensions
imgdim = abs(diff([bb [1;1]]*inv(mat')));
imgdim = round(imgdim(1:3))+1;


%%% Create the transformation matrix from the input image to the desired
%%% output image.
tm = mat\h1(1).mat;

%%% Use spm_slice_vol to reslice the image to the desired specifications
NN = zeros([imgdim numel(h1)]);
% persisText;
for jj = 1:length(h1);
    
    V = h1(jj);
    V = rmfield(V,'private');
    V.fname = strtrim([prefix V.fname]);
    V.dim = imgdim(1:3);
    V.mat = mat;
%     V = spm_create_vol(V);
     
%     persisText(['Processing Volume #' num2str(jj)]);
    N = zeros(imgdim);
    for i = 1:imgdim(3)
        M = tm\spm_matrix([0 0 i]);
        img = spm_slice_vol(h1(jj), M, imgdim(1:2), interpOrd);
        N(:,:,i) = img;
%         if ismask
%             img = round(img);
%         end
%         spm_write_plane(V, img, i);
    end
    NN(:,:,:,jj) = N;
    if ~isempty(prefix)
        spm_write_vol(V,N);
    end
end

if numel(NN) > numel(N)
    N = NN;
end

% persisText
% disp('Done.')

%%% And we are done!

