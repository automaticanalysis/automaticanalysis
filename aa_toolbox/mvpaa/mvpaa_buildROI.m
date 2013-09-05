function [ indROI voxels ] = mvpaa_buildROI( centreROI, subROI, brainSize)
%BUILD_ROI Builds an ROI mask
% - The ROI sphere is centred on one central voxel each.
% - The subROI is by default an array with the base indices of an ROI
% centred at [0 0 0], making many of them negative
% - The brainSize is by default a vector with the values for hte 3 dims

if length(centreROI) ~= 3
    error('Your ROI centre are not correctly specified')
end

include = 0;
% Adapt indices to the centres we wish to take
for a = 1:3 % 3 dims
    subROI(:,a) = subROI(:,a) + centreROI(a);
    % Check which are less than 1 or more than their respective brain size
    include = include + (subROI(:,a) > 0 & subROI(:,a) <= brainSize(a));
end

include = include == length(centreROI);
voxels = sum(include);

% Remove voxels outside the matrix
TsubROI = nan(voxels,3);
for a = 1:3 % 3 dims
    TsubROI(:,a) = subROI(include,a);
end
       
% Get the indices...
indROI = sub2ind(brainSize, TsubROI(:,1), TsubROI(:,2), TsubROI(:,3));   