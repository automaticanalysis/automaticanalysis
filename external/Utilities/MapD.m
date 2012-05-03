function MapD
% For use with two-sample t-tests
load I;
i1 = find(I.X(:,1)==1);
i2 = find(I.X(:,2)==1);

[m1 h] = openIMG(char(I.Scans(i1)));
m1 = reshapeWholeBrain(size(m1),m1);

m2 = openIMG(char(I.Scans(i2)));
m2 = reshapeWholeBrain(size(m2),m2);

D = cohensD(m1,m2);

hh = h(1);
hh.fname = 'CohensD.nii';
V = zeros(hh.dim);
V(:) = D;

spm_write_vol(hh,V);
