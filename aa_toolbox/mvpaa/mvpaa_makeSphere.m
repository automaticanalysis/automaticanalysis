function [ sphere voxels ] = mvpaa_makeSphere( radius )
%MAKE_SPHERE Makes a sphere
% The sphere is centred on one central voxel.
% The radius excludes the central voxel, so diameter is radius*2 + 1.

sphere=zeros(radius*2 + 1, radius*2 + 1,radius*2 + 1);
[X Y Z] = meshgrid(1:(radius*2 + 1), 1:(radius*2 + 1), 1:(radius*2 + 1));
X = X - (radius + 1);
Y = Y - (radius + 1);
Z = Z - (radius + 1);
D = sqrt(X.^2 + Y.^2 + Z.^2);

sphere(D < (radius + 0.001)) = 1;
voxels = sum(sum(sum(sphere)));