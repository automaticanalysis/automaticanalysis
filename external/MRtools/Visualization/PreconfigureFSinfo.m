%%
%%% Written by Aaron P. Schultz - aschultz@martinos.org
%%%
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
%% Create FSaverage mat file from freesurfer source
% Num Verts
% fsaverage3 = 642
% fsaverage4 = 2562
% fsaverage5 = 10242
% fsaverage6 = 40962
% fsaverage  = 163842

clear all
pth = ['/autofs/space/plato_002/users/freesurfer/subjects/fsaverage/'];
[lVert, lFace] = read_surf([pth 'surf/lh.inflated_avg']);
[rVert, rFace] = read_surf([pth 'surf/rh.inflated_avg']);
T.inflated.lVert = lVert;
T.inflated.lFace = lFace+1;
T.inflated.rVert = rVert;
T.inflated.rFace = rFace+1;

[lVert, lFace] = read_surf([pth 'surf/lh.pial_avg']);
[rVert, rFace] = read_surf([pth 'surf/rh.pial_avg']);
T.pial.lVert = lVert;
T.pial.lFace = lFace+1;
T.pial.rVert = rVert;
T.pial.rFace = rFace+1;

[lVert, lFace] = read_surf([pth 'surf/lh.white_avg']);
[rVert, rFace] = read_surf([pth 'surf/rh.white_avg']);
T.white.lVert = lVert;
T.white.lFace = lFace+1;
T.white.rVert = rVert;
T.white.rFace = rFace+1;

[lSulc, fnum] = read_curv([pth 'surf/lh.avg_sulc']);
[lCurv, fnum] = read_curv([pth 'surf/lh.avg_curv']);
[rSulc, fnum] = read_curv([pth 'surf/rh.avg_sulc']);
[rCurv, fnum] = read_curv([pth 'surf/rh.avg_curv']);
[rThk, fnum] = read_curv([pth 'surf/lh.avg_thickness']);
[lThk, fnum] = read_curv([pth 'surf/rh.avg_thickness']);

T.lSulc = lSulc;
T.lCurv = lCurv;
T.lThk = lThk;
T.rSulc = rSulc;
T.rCurv = rCurv;
T.rThk = rThk;

lh_map = dlmread([pth 'label/lh.cortex.label'],' ',2,0);
rh_map = dlmread([pth 'label/rh.cortex.label'],' ',2,0);

lMNI = lh_map(:,3:2:7);
lv = lh_map(:,1);

rMNI = rh_map(:,3:2:7);
rv = rh_map(:,1);

T.map.lMNI = lMNI;
T.map.rMNI = rMNI;
T.map.lv = lv;
T.map.rv = rv;

save fsaverage.mat T

%% Create a mapping file for prettier surface images, this creates a wieghted average mapping from the volume specified to the fsaverage based on voxel proximity to vertices.
%%% Code from the first section is not neceassry for making the mapping
%%% file if you already have the fsaverage mat files.

%%% This can take a little while (at least a few minutes), but need only be
%%% done once per image size.  (Note size is an orientation, with a given
%%% voxel size, with a given bounding box.);
clear
%%% Choose a surface to map to
load fsaverage6.mat
%%% Choose an image file to map
input  = '/autofs/space/schopenhauer_003/users/BigTemplate/Orthomax_NormOn_20/Stuff/DMN.nii';
[m h] = openIMG(input);
%%% Choose a proximity threshold in mm
prox = 4.5;

[x y z] = ind2sub(h.dim,(1:numel(m))');
mat = [x y z ones(numel(z),1)];
mni = mat*h.mat';
mni = mni(:,1:3);

he = h;

nn = 100;
lWeights = zeros(size(T.map.lMNI,1),nn);
lVoxels = ones(size(T.map.lMNI,1),nn);
for ii = 1:size(T.map.lMNI,1);
    dd = sqrt(sum((repmat(T.map.lMNI(ii,:),size(mni,1),1)-mni).^2,2));
    
    i1 = find(dd<prox);
    lVoxels(ii,1:numel(i1))=i1;
    tmp = prox-dd(i1);
    lWeights(ii,1:numel(i1))=tmp./sum(tmp);
end


rWeights = zeros(size(T.map.rMNI,1),nn);
rVoxels = ones(size(T.map.rMNI,1),nn);
for ii = 1:size(T.map.rMNI,1);
    dd = sqrt(sum((repmat(T.map.rMNI(ii,:),size(mni,1),1)-mni).^2,2));
    
    i1 = find(dd<prox);
    rVoxels(ii,1:numel(i1))=i1;
    tmp = prox-dd(i1);
    rWeights(ii,1:numel(i1))=tmp./sum(tmp);
end

i1 = find(sum(lWeights)~=0);
lWeights = lWeights(:,i1);
lVoxels = lVoxels(:,i1);

i1 = find(sum(rWeights)~=0);
rWeights = rWeights(:,i1);
rVoxels = rVoxels(:,i1);

lVoxels(lVoxels==0)=1;
rVoxels(rVoxels==0)=1;

MP.lVoxels = lVoxels;
MP.lWeights = lWeights;
MP.rVoxels = rVoxels;
MP.rWeights = rWeights;
MP.header = h;
header = h;

%%% save the mapping info to a mat file that can be use as an input ro
%%% surfPlot.m

save Test.mat MP header



    
