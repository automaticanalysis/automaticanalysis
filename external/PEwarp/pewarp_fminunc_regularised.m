function x = pewarp_fminunc_regularised(targetFile, sourceFile, order, xscl) % maskFile

global bX bY dbY bZ target8 source mask krn xscale regstrength

% histogram smoothing
fwhm = 4;

regstrength = single(.1);

if nargin < 3
    order = [10 16 10];
end

if nargin < 4
    xscale = single(1e6);
else
    xscale = single(xscl);
end

targetVol = spm_vol(targetFile);

target = zeros(targetVol.dim);
for i = 1 : targetVol.dim(3)
   target(:, :, i) = spm_slice_vol(targetVol, spm_matrix([0 0 i]), targetVol.dim(1 : 2), 1);
end;

targetMin = min(target(:));
targetMax = max(target(:));

target8 = uint8(round((target - targetMin) * 255 / (targetMax - targetMin)));

sourceVol = spm_vol(sourceFile);

source = single(zeros(sourceVol.dim));
for i = 1 : sourceVol.dim(3);
    source(:, :, i) = spm_slice_vol(sourceVol, spm_matrix([0 0 i]), sourceVol.dim(1 : 2), 1);
end;

sourceMin = min(source(:));
sourceMax = max(source(:));

source = (source - sourceMin) * 255 / (sourceMax - sourceMin);

%maskVol = spm_vol(maskFile);

%mask = uint8(zeros(maskVol.dim));
%for i = 1 : maskVol.dim(3);
%    mask(:, :, i) = spm_slice_vol(maskVol, spm_matrix([0 0 i]), maskVol.dim(1 : 2), 1);
%end;

mask = uint8(ones(size(source)));

bX = single(spm_dctmtx(sourceVol.dim(1), order(1)));
bY = single(spm_dctmtx(sourceVol.dim(2), order(2)));
dbY = single(spm_dctmtx(sourceVol.dim(2), order(2), 'diff'));
bZ = single(spm_dctmtx(sourceVol.dim(3), order(3)));

% From spm_coreg
% lim  = ceil(2 * fwhm);
lim = ceil(fwhm);
krn = smoothing_kernel(fwhm, -lim : lim) ;
krn = single(krn / sum(krn));

x0 = zeros(prod(order), 1);

% x = fminsearch(@pewarp_fminunc_costwrapper, x0, optimset('Display', 'iter-detailed',  'TolFun', 1e-10));
x = fminunc(@pewarp_fminunc_regularised_costwrapper, x0, optimset('Display', 'iter-detailed', ...
    'OutputFcn', @pewarp_fminunc_regularised_plot, 'LargeScale', 'off', 'Diagnostics', 'on', ...
    'FinDiffType', 'central', 'MaxFunEvals', 1000000, 'TolX', 1e-10));
% x = pewarp_spm_powell(x, xi, tolsc, maxiter, 'mex_pewarpcost', bX, bY, dbY, bZ, target8, source, mask, krn, xscale);


% From spm_coreg
function krn = smoothing_kernel(fwhm,x)

% Variance from FWHM
s = (fwhm/sqrt(8*log(2)))^2+eps;

% The simple way to do it. Not good for small FWHM
% krn = (1/sqrt(2*pi*s))*exp(-(x.^2)/(2*s));

% For smoothing images, one should really convolve a Gaussian
% with a sinc function.  For smoothing histograms, the
% kernel should be a Gaussian convolved with the histogram
% basis function used. This function returns a Gaussian
% convolved with a triangular (1st degree B-spline) basis
% function.

% Gaussian convolved with 0th degree B-spline
% int(exp(-((x+t))^2/(2*s))/sqrt(2*pi*s),t= -0.5..0.5)
% w1  = 1/sqrt(2*s);
% krn = 0.5*(erf(w1*(x+0.5))-erf(w1*(x-0.5)));

% Gaussian convolved with 1st degree B-spline
%  int((1-t)*exp(-((x+t))^2/(2*s))/sqrt(2*pi*s),t= 0..1)
% +int((t+1)*exp(-((x+t))^2/(2*s))/sqrt(2*pi*s),t=-1..0)
w1  =  0.5*sqrt(2/s);
w2  = -0.5/s;
w3  = sqrt(s/2/pi);
krn = 0.5*(erf(w1*(x+1)).*(x+1) + erf(w1*(x-1)).*(x-1) - 2*erf(w1*x   ).* x)...
      +w3*(exp(w2*(x+1).^2)     + exp(w2*(x-1).^2)     - 2*exp(w2*x.^2));

krn(krn<0) = 0;
return;


