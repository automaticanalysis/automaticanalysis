function [N h] = resizeVol3(h,targ,ord)

if ischar(h)
    h = spm_vol(h);
end

if ischar(targ)
    targ = spm_vol(targ);
end

if nargin<3 || isempty(ord)
    ord = [3 0];
end

if numel(ord)<2
    ord(2) = 0;
end
FillVal=ord(2);

m = spm_read_vols(h);

d = [ord(1) ord(1) ord(1) 0 0 0];
c = spm_bsplinc(m,d);


[x1,x2,x3] = ndgrid(1:targ.dim(1),1:targ.dim(2),1:targ.dim(3));
[msk,y1,y2,y3] = getmask(inv(targ.mat\h.mat),x1,x2,x3,targ.dim(1:3),[0 0 0]);
N = spm_bsplins(c,y1,y2,y3,d);
if FillVal~=0
    N(msk==0)=FillVal;
end

% targ.fname = 'reori.nii';
% targ.dt = h.dt;
% spm_write_vol(targ,N);
%%
% X1 = ones(numel(m1),4);
% [X1(:,1), X1(:,2), X1(:,3)] = ind2sub(size(m1),1:numel(m1));
% oX1=X1;
% X1 = (X1*R');
%% Straight resampling is faster
% clc
% [m h] = openIMG('02_DefaultMode.nii'); 
% 
% tic
% [oX,oY,oZ] = meshgrid(1:size(m,2),1:size(m,1),1:size(m,3));
% [nX,nY,nZ] = meshgrid(1:.5:size(m,2),1:.5:size(m,1),1:.5:size(m,3));
% Vq = interp3(oX,oY,oZ,m,nX,nY,nZ,'cubic',0);
% toc

%% try also using permute to for orthongal roations.
% [m1 h1] = openIMG('defaultUnderlay.nii');
% h = h1;
% % ord = [1 3 2];
% % ord = [2 1 3];
% % ord = [2 3 1];
% % ord = [3 1 2];
% ord = [3 2 1];
% m2 = permute(m1,ord);
% 
% % h1.mat(:,1) = h.mat(:,ord(1));
% % h1.mat(:,2) = h.mat(:,ord(2));
% % h1.mat(:,3) = h.mat(:,ord(3));
% h1.mat = h1.mat(:,[ord 4]);
% 
% 
% h.mat\h1.mat
% 
% return
% h1.fname = 'A1.nii';
% h1.dim = h1.dim(ord);
% spm_write_vol(h1,m2);
% 
% 
% FIVE('A1.nii');
% 
% 
% %%
% figure(99); clf; subplot(1,2,1); imagesc(squeeze(m1(:,:,50))); subplot(1,2,2); imagesc(squeeze(m2(:,:,50))); shg
% %%
% figure(99); clf; subplot(1,2,1); imagesc(squeeze(m1(:,50,:))); subplot(1,2,2); imagesc(squeeze(m2(50,:,:))); shg
% %%
% figure(99); clf; subplot(1,2,1); imagesc(squeeze(m1(50,:,:))); subplot(1,2,2); imagesc(squeeze(m2(:,50,:))); shg
% 


