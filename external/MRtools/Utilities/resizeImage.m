function [mm h] = resizeImage(input,voxdim,filename)
% voxdim = [1 2 3];

if ischar(input)
    [m h] = openIMG(input); 
elseif isstruct(input)
    if isfield(input,'fname')
        h = input;
        m = spm_read_vols(input);
    else
        h = input.h;
        m = input.m;
    end
end

a = spm_imatrix(h.mat);
adj = voxdim./abs(a(7:9));

lens = ceil(size(m)./adj);

nx = linspace(1,size(m,2),lens(2));
ny = linspace(1,size(m,1),lens(1));
nz = linspace(1,size(m,3),lens(3));


[oX,oY,oZ] = meshgrid(1:size(m,2),1:size(m,1),1:size(m,3));
[nX,nY,nZ] = meshgrid(nx,ny,nz);
mm = interp3(oX,oY,oZ,m,nX,nY,nZ,'cubic',0);

a(7:9)=voxdim.*sign(a(7:9));
h.mat = spm_matrix(a);
h.dim = lens;
h.fname = [];

if nargin>2 && ~isempty(filename)
    h.fname = filename;
    spm_write_vol(h,mm)
end