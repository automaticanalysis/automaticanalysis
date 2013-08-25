function H = spm_hist2_weighted(reference,source,M,sg,wgt)
% 
% H = spm_hist2_weighted(VG.uint8,VF.uint8, VF.mat\spm_matrix(x(:)')*VG.mat ,sg, wgt)
%
% A custom version of spm_hist2 to be able to perform weighted
% coregistration.
%
% spm_hist was originally called by:
% H = spm_hist2(VG.uint8,VF.uint8, VF.mat\spm_matrix(x(:)')*VG.mat ,sg);
% 
% this one is the same only a weighting volume can be given. it should have
% the same orientation and dimensions as the reference volume
%


% the approach of spm_hist2 is simulated, but is adapted at some points.
H = zeros(256);

dims = size(reference);

[Y,X,Z] = meshgrid(1:sg(2):dims(2),1:sg(1):dims(1),1:sg(3):dims(3));    % sample data

X = X + rand(size(X))*sg(1);    % jitter data
Y = Y + rand(size(Y))*sg(2);
Z = Z + rand(size(Z))*sg(3);

% x coordinates of the source dataset (which is rotated/translated)
X_new = M(1,1)*X + M(1,2)*Y + M(1,3)*Z + M(1,4);      
Y_new = M(2,1)*X + M(2,2)*Y + M(2,3)*Z + M(2,4);
Z_new = M(3,1)*X + M(3,2)*Y + M(3,3)*Z + M(3,4);

vg = spm_sample_vol(reference,X,Y,Z,1);
ivg = floor(vg+0.5);    % faster than round.m

vf = spm_sample_vol(source,X_new,Y_new,Z_new,1);
ivf = floor(vf);
diff_vf = vf - ivf;

wgt = spm_sample_vol(wgt,X,Y,Z,1);

% first part
for cnt = 1:size(vg,1)*size(vg,2);
    H(ivf(cnt)+1,ivg(cnt)+1) = H(ivf(cnt)+1,ivg(cnt)+1) + (1 - diff_vf(cnt))*wgt(cnt);
end

% second part in a separate for loop to avoid an if statement within a for
% loop
ind = find(ivf<255);
for cnt = 1:length(ind)
    H(ivf(ind(cnt))+2,ivg(ind(cnt))+1) = H(ivf(ind(cnt))+2,ivg(ind(cnt))+1) + diff_vf(ind(cnt))*wgt(ind(cnt));
end



