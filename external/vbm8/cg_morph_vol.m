function vol = cg_morph_vol(in,action,n,th,vx_vol);
% morphological operations to 3D data
%__________________________________________________________________________
% Christian Gaser
% $Id: cg_morph_vol.m 404 2011-04-11 10:03:40Z gaser $

rev = '$Rev: 404 $';

if nargin < 5, vx_vol = [1 1 1]; end
if nargin < 4, th = 0.5; end
if nargin < 3, n = 1; end
if nargin < 2, action = 'open'; end
if nargin < 1
	error('No arguments given.');
end

th = th*max(double(in(:)));

vol = uint8(in>th);

kx = [1 1 1];
ky = [1 1 1];
kz = [1 1 1];

order = sum(kx(:) ~= 0)*sum(ky(:) ~= 0);


switch lower(action)
	case 'dilate'
	%=======================================================================
	% enlarge image according to number of dilations
	sz = size(vol);
	vol2 = zeros(sz(1)+(2*n),sz(2)+(2*n),sz(3)+(2*n),'uint8');
	vol2(n+1:sz(1)+n,n+1:sz(2)+n,n+1:sz(3)+n) = vol;
	for i = 1:n
		spm_conv_vol(vol2,vol2,kx,ky,kz,-[1 1 1]);
		vol2 = uint8(vol2~=0);
	end
	vol = vol2(n+1:sz(1)+n,n+1:sz(2)+n,n+1:sz(3)+n);
	clear vol2

	case 'erode'
	%=======================================================================
	for i = 1:n
		spm_conv_vol(vol,vol,kx,ky,kz,-[1 1 1]);
		vol = uint8(vol>=order);
	end

	case 'close'
	%=======================================================================
	% enlarge image according to number of dilations
	sz = size(vol);
	vol2 = zeros(sz(1)+(2*n),sz(2)+(2*n),sz(3)+(2*n),'uint8');
	vol2(n+1:sz(1)+n,n+1:sz(2)+n,n+1:sz(3)+n) = vol;
	for i = 1:n
		spm_conv_vol(vol2,vol2,kx,ky,kz,-[1 1 1]);
		vol2 = uint8(vol2~=0);
	end
	for i = 1:n
		spm_conv_vol(vol2,vol2,kx,ky,kz,-[1 1 1]);
		vol2 = uint8(vol2>=order);
	end
	vol = vol2(n+1:sz(1)+n,n+1:sz(2)+n,n+1:sz(3)+n);
	clear vol2
	
	case 'open'
	%=======================================================================
	for i = 1:n
		spm_conv_vol(vol,vol,kx,ky,kz,-[1 1 1]);
		vol = uint8(vol>=order);
	end
	for i = 1:n
		spm_conv_vol(vol,vol,kx,ky,kz,-[1 1 1]);
		vol = uint8(vol~=0);
  end
  
  case 'labclose'
	%=======================================================================
  [ROI,num] = spm_bwlabel(double(~vol),6);
  num       = hist( ROI( ROI(:)>0 ) , 1:num);
  [tmp,num]   = max(num(1:end));  
  vol = uint8(1 - (ROI==num));	
  
  case 'labopen'
	%=======================================================================
  [ROI,num] = spm_bwlabel(double(vol),6);
  num       = hist( ROI( ROI(:)>0 ) , 1:num);
  [tmp,num]   = max(num(1:end));  
  vol = uint8(ROI==num);
  
  case 'distclose'
  %=======================================================================
  vol = vbdist(single(vol),true(size(vol)),vx_vol);
  vol = vbdist(single(vol>n),true(size(vol)),vx_vol)>n;
  
  case 'distopen'
  %=======================================================================
  vol = vbdist(1-single(vol),true(size(vol)),vx_vol);
  vol = vbdist(single(vol>n),true(size(vol)),vx_vol)<n;
  
	otherwise
		error('Unknown action');
end

if isa(in,'double')
	vol = double(vol);
end

if isa(in,'single')
	vol = single(vol);
end

return;
