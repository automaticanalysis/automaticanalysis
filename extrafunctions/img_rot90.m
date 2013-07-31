% Rotate a 3D image along the axis perpenidcular to the x-y plane with 90 degrees.
% Tibor Auer MRC CBU Cambridge 2012-2013

function fo = img_rot90(fi,n)
if nargin < 2, n = 1; end
nslice = size(fi,3);
for i = 1:nslice
    fo(:,:,i) = rot90(fi(:,:,i),n);
end
