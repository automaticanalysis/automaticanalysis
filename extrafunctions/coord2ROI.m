%% coord2ROI
% This function creates an ROI nifti file, using as inputs 4 variables:
% 1) coord    -- coordinate in the template space (e.g. [10, 15, 17] mm in MNI)
% 2) radius   -- radius in the template space (e.g. 5 mm in MNI)
% 3) ROIname  -- name of the ROI (with '.nii' or '.img' extension)
% 4) template -- filename of the MNI or Talaraich template to use

function coord2ROI(coord, radius, ROIname, template)

if nargin < 1
    coord = [0 0 0]; % 46.0 64.0 37.0
end
if nargin < 2
    radius = 5; % in millimetres
end
if nargin < 3
    ROIname = 'myROI.nii';
end
if nargin < 4
    template = fullfile(spm('dir'), 'templates/T1.nii');
end

% Put coordinate vector in right orientation...
coord = coord(:);

% Get the template...
V = spm_vol(template);

% Transform coordinate system into voxels...
coord = V.mat \ [coord; 1];
coord = coord(1:3);

% Get the number of mm per voxel...
mmVox = vox2mm(V);

% Make radius in voxel space...
radius = [radius radius radius];
radius = radius ./ mmVox;

% Create a meshgrid showing which voxels are here...
[Y X Z] = meshgrid(1:V.dim(2), 1:V.dim(1), 1:V.dim(3));

% Subtract from meshgrids our coordinates...
X = X - coord(1);
Y = Y - coord(2);
Z = Z - coord(3);

% ...and divide by our radius
X = X ./ radius(1);
Y = Y ./ radius(2);
Z = Z ./ radius(3);

% Get absolute distance
M = sqrt(X.^2 + Y.^2 + Z.^2);

% ...and mask
M = M <= 1;

% Write the mask
V.dt = [2 0];
V.fname = ROIname;

spm_write_vol(V,M);