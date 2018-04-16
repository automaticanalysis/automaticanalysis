function k = smoothing_kernel(s,VOX)
% P = spm_vol(char(P));
% VOX = sqrt(sum(P(1).mat(1:3,1:3).^2));

s  = s./VOX;                        % voxel anisotropy
s1 = s/sqrt(8*log(2));              % FWHM -> Gaussian parameter

x  = round(6*s1(1)); x = -x:x; x = spm_smoothkern(s(1),x,1); x  = x/sum(x);
y  = round(6*s1(2)); y = -y:y; y = spm_smoothkern(s(2),y,1); y  = y/sum(y);
z  = round(6*s1(3)); z = -z:z; z = spm_smoothkern(s(3),z,1); z  = z/sum(z);

i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;


k = ones(numel(x),numel(y),numel(z));
[i1 i2 i3] = ind2sub(size(k),1:numel(k));
val = prod([x(i1)' y(i2)' z(i3)'],2);
k(:) = val;