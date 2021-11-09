function varargout = read_voxels_matlab(varargin)

%   read_voxels: 
%       inputs: header, coordinates (n x 3 (XYZ))
%       output: 1 x n vector of voxel values

coord = num2cell(varargin{2});

img = read_image_matlab(varargin{[1 3:end]});

varargout{1} = arrayfun(@(vox) img(coord{vox,:}), 1:size(coord,1));