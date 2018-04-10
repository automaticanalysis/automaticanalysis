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

%%% How to use surfPlot outside of the FIVE GUI
clear all; %close all; clc

%%% Specify an input image
obj.input = '/autofs/space/schopenhauer_003/users/BigTemplate/Orthomax_NormOn_20/Stuff/DMN.nii';
%%% or specify a header plus 3D matrix
% [m h] = openIMG('/autofs/space/schopenhauer_003/users/BigTemplate/Orthomax_NormOn_20/Stuff/DMN.nii');
% obj.input.m = m;
% obj.input.he = h;


obj.figno = 21; % Figure number for output plot
obj.newfig = 1; % Specifies whether or not to re-render the surface or to add to existing render
obj.overlaythresh = [70];  % Specifes the threshold to be applies
% obj.overlaythresh = [-70 70];  % Can also be specified as two sided threshold
obj.colorlims = [0 inf];  % Set the color limits to be used inf and -inf are translated to min and max values
obj.colomap = 'jet';  % Choose your color map:  See colmap.m for options
obj.direction = '+';  % sets direction of one sided threshold + = greater than -= less than
obj.reverse = 0;  % Option to reverse the image; i.e. m*-1

%%% Use either these:
% obj.mappingfile = []; % if no mapping file is specifed mapping is done by an average of the nearest neighbors
% obj.nearestneighbor = 0; % if = 1, only the value from the closest voxel will be used, useful for maskings and label images
%%% Or specify a mapping file.
obj.mappingfile = 'Test.mat';  %%% See PreconfigureFSinfo.m for an example of how to create a mapping file.

obj.round = 0;  % if = 1, rounds all values on the surface to nearest whole number.  Useful for masks

obj.fsaverage = 'fsaverage6';  %% Set which fsaverage to map to e.g. fsaverage, fsaverage3, fsaverage6
obj.surface = 'pial';          %% Set the surface: inflated, pial, or white
obj.shading = 'curv';          %% Set the shading information for the surface: curv, sulc, or thk
obj.shadingrange = [.1 .7];    %% Set the min anx max greyscale values for the surface underlay (range of 0 to 1)
obj.Nsurfs = -1;              %% Choose which hemispheres and surfaces to show:  4=L/R med/lat;  2= L/R lat; 1.9=L med/lat; 2.1 = R med/lat; -1= L lat; 1-R lat;

%%% h1 =  underlay handles
%%% hh1 = overlay handles
[h1 hh1] = surfPlot(obj);

% set(gcf,'color','w'); shg

%% Add a second overlay to the first endering
%%% A second overlay can then be added by changing a few fields of obj and
%%% passing obj to surfPlot.
obj.input = '/autofs/space/schopenhauer_003/users/BigTemplate/Orthomax_NormOn_20/Stuff/FPCN.nii';
obj.colomap = 'cool';
obj.newfig=0;
obj.overlaythresh = [50];
obj.colorlims = [50 inf];
[h2 hh2] = surfPlot(obj);
shg

% set(hh1,'FaceAlpha',.85)
% set(hh2,'FaceAlpha',.85)
