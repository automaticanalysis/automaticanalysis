function [g,w,c] = cg_cleanup_gwc(g,w,c, level)
% use morphological operations to cleanup GM/WM/CSF
%__________________________________________________________________________
% Christian Gaser
% $Id: cg_cleanup_gwc.m 404 2011-04-11 10:03:40Z gaser $

rev = '$Rev: 404 $';

if nargin<4, level = 1; end;

th0 = 0.6;

mx_w = max(double(w(:)));

if level==2
	th1 = 0.2*mx_w/255;
else
	th1 = 0.15*mx_w/255;
end;

th0 = th0*mx_w/255;

% use only largest WM cluster at th0 as seed parameter
% mask out all voxels where wm is lower than gm or csf
wmask = w;
wmask(find((double(w) < double(c)) & (double(w) < double(g)))) = 0;
b = mask_largest_cluster(wmask, th0);

% Build a 3x3x3 seperable smoothing kernel
%-----------------------------------------------------------------------
kx=[0.75 1 0.75];
ky=[0.75 1 0.75];
kz=[0.75 1 0.75];
sm=sum(kron(kron(kz,ky),kx))^(1/3);
kx=kx/sm; ky=ky/sm; kz=kz/sm;

% Erosions and conditional dilations
%-----------------------------------------------------------------------

niter = 45;
spm_progress_bar('Init',niter,'Clean up','Iterations completed');

for j=1:niter,
        if j>2, th=th1; else th=th0; end; % Dilate after two its of erosion.
        for i=1:size(b,3),
                gp = double(g(:,:,i));
                wp = double(w(:,:,i));
                bp = double(b(:,:,i))/255;
                bp = (bp>th).*(wp+gp);
                b(:,:,i) = uint8(round(bp));
        end;
        spm_conv_vol(b,b,kx,ky,kz,-[1 1 1]);
        spm_progress_bar('Set',j);
end;

th = 0.05;

b2 = double(w);
b2 = b2 + double(g);
b2 = b2 + double(c);
b2 = cg_morph_vol(b2,'open',1,th1);
b2 = mask_largest_cluster(b2,0.5);
b2 = double(cg_morph_vol(b2,'close',2,0.5));

for i=1:size(b,3),
        gp       = double(g(:,:,i))/255;
        wp       = double(w(:,:,i))/255;
        cp       = double(c(:,:,i))/255;
        bp       = double(b(:,:,i))/255;
        bp       = ((bp>th).*(wp+gp))>th;
        gwc      = gp + wp + cp + eps;
        g(:,:,i) = uint8(round(255*gp.*bp./gwc));
        w(:,:,i) = uint8(round(255*wp.*bp./gwc));
        c(:,:,i) = uint8(round(255*cp.*b2(:,:,i)./gwc));
end;
return

%=======================================================================
function y = mask_largest_cluster(y, th)

if nargin < 2
	th = 0.5;
end

sz = size(y);

th = th*max(double(y(:)));

mask = y > th;
Q = find(mask);

Qth = find(y <= th & y>0);
yth = y(Qth);

% save memory by using bounding box where y > th
[indx, indy, indz] = ind2sub(size(mask),Q);
indx = max((min(indx) - 1),1):min((max(indx) + 1),sz(1));
indy = max((min(indy) - 1),1):min((max(indy) + 1),sz(2));
indz = max((min(indz) - 1),1):min((max(indz) + 1),sz(3));

[A,num] = spm_bwlabel(double(mask(indx,indy,indz)),26);

clear mask

% interrupt if cluster was > 7.5% of whole image to save time
max_A = max(A(:));
sz_cluster = zeros(max_A,1);
for i=1:max_A
	QA = find(A == i);
	ind = i;
	if length(QA)/prod(size(A)) > 0.075
		break
	end
	sz_cluster(i) = length(QA);
end

if length(QA)/prod(size(A)) <= 0.075
	[mx, ind] = max(sz_cluster);
	QA = find(A == ind);
end

QA0 = find(A ~= ind);
A = y(indx,indy,indz);
A(QA0) = 0;

y(indx,indy,indz) = A;
y(Qth) = yth;


spm_progress_bar('Clear');
return;
