function varargout = read_voxels_matlab(varargin)

%   read_voxels: 
%       inputs: header, coordinates (n x 3 (XYZ))
%       output: 1 x n vector of voxel values

coord = num2cell(varargin{2});
if size(coord,1) == 1 % concatenated in a wrong way
    nVox = size(coord,2)/3;
    coord = reshape(coord,nVox,3);
end

img = read_image_matlab(varargin{[1 3:end]});

varargout{1} = arrayfun(@(vox) img(coord{vox,:}), 1:size(coord,1));