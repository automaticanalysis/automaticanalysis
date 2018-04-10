function [nm nh] = canonOri(m,h,refH)
%%% m is the data matrix to be rotated/
%%% h is the nifti style header for m.
%%% refH is an optional input of a second header to use to define the image
%%% size.
%%%
%%% nm is the output rotated/resliced matrix
%%% nh is configured header for nm.
%%%

orig = spm_imatrix(h.mat);

to = spm_matrix([0 0 0 orig(4:6) sign(orig(7:9)) 0 0 0]);
rm = (diag([-1 1 1])\to(1:3,1:3));
[tr order] = sort([1 2 3] * rm');
voxdim = abs(orig(6+order));

%%%

if nargin>2
    if isstruct(refH)
        bb = world_bb(refH);
    else
        bb = refH;
    end
else
    bb = world_bb(h);
end

TM = diag([voxdim .* [-1 1 1] 1]);
TM(1,4) = bb(2,1)+voxdim(1);
TM(2,4) = bb(1,2)-voxdim(2);
TM(3,4) = bb(1,3)-voxdim(3);

ms = round(([bb(2,:) 1]*inv(TM'))+([bb(1,:) 1]*inv(TM')))-1;

nh = h; 
nh.fname = [];
nh.mat = TM; 
nh.dim = ms(1:3); 
%%%
rm(abs(rm)<1e-6)=0;
if all(sum(abs(rm),1)==1) && all(sum(abs(rm),2)==1) && nargin<3
    
    [r c] = find(abs(rm') == 1);
    if all(r==[1 2 3]')
        nm = m;
        disp('Data is good to go'); %return
    else
        disp('Shuffling Dimensions');
        nm = permute(m,r');
    end
    
    for ii = 1:size(rm,2)
        if sum(rm(ii,:))<0
            disp('Flipping Dimensions');
            nm = flip(nm,ii);
        end
    end
else
    disp('ReslicingData'); % return
    nm = resizeVol4(m,h,nh,[3 0]);
end
%%%
% nh.fname = 'ATest.nii';
% spm_write_vol(nh,nm);
% FIVE('ATest.nii');
%%
% figure(99); clf;
% mm = openIMG('defaultUnderlay.nii');
% 
% i1 = round(size(m)/2);
% subplot(3,3,1); imagesc(squeeze(m(i1(1),:,:)));  axis equal; shg
% subplot(3,3,2); imagesc(squeeze(m(:,i1(2),:))); axis equal; shg
% subplot(3,3,3); imagesc(squeeze(m(:,:,i1(3))));  axis equal; shg
% 
% i1 = round(size(mm)/2);
% subplot(3,3,4); imagesc(squeeze(mm(i1(1),:,:)));  axis equal; shg
% subplot(3,3,5); imagesc(squeeze(mm(:,i1(2),:)));  axis equal; shg
% subplot(3,3,6); imagesc(squeeze(mm(:,:,i1(3))));   axis equal; shg
% 
% i1 = round(size(nm)/2);
% subplot(3,3,7); imagesc(squeeze(nm(i1(1),:,:)));  axis equal; shg
% subplot(3,3,8); imagesc(squeeze(nm(:,i1(2),:)));  axis equal; shg
% subplot(3,3,9); imagesc(squeeze(nm(:,:,i1(3))));  axis equal; shg
%%
