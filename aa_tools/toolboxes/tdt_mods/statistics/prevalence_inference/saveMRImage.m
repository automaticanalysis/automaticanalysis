function V = saveMRImage(data, filename, mat, descrip)

% saves a 3d-matrix to an MR image file, using SPM routines
%
% V = saveMRImage(data, filename, mat = eye(4), descrip = 'saveMRImage')
%
% data:         image data to be written, 3d matrix
% filename:     name of the file to write to
% mat:          transformation from voxel indices to mm coordinates (aka v2mm)
% descrip:      string describing contents
% V:            spm volume struct
%
% The file format is determined by the filename extension.
% The data type is is determined by the class of data.

% Author: probably Carsten Allefeld, 2014/4/7

if nargin < 3
    mat = eye(4);
end

if nargin < 4
    descrip = 'saveMRImage';
end

if size(data, 4) > 1
    warning('only the first volume is written')
end

if isfloat(data) && ~isreal(data)
    warning('only real part of complex data is written')
    data = real(data);
end

% determine data type
dt = class(data);
if isa(data, 'single')
    dt = 'float32';
elseif isa(data, 'double')
    dt = 'float64';
end
dt = spm_type(dt);
if isnan(dt)
    error('data type is not supported')
end

V = struct;
V.fname = filename;
V.dim = size(data,1:3);
V.dt = [dt 0];
V.mat = mat;
V.descrip = descrip;

V = spm_write_vol(V, data);

if nargout == 0
    clear V
end
