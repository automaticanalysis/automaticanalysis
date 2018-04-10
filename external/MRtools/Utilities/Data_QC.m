function S = Data_QC(fn,thresh,ROIs)
%%% This script replicates stackcheck_nifti.  It will produce Mean, SD, and
%%% SNR images as well as computing and writing out Global mean, SD, and SNR
%%% values.
%%%
%%% fn:     A nifti file name.
%%%
%%% Written by Aaron Schultz (aschultz@martinos.org)
%%%
%%% Last Updated Dec. 11 2012;
%%%
%%% Copyright (C) 2012,  Aaron P. Schultz
%%%
%%% Supported in part by the NIH funded Harvard Aging Brain Study (P01AG036694) and NIH R01-AG027435 
%%%
%%% This program is free software: you can redistribute it and/or modify
%%% it under the terms of the GNU General Public License as published by
%%% the Free Software Foundation, either version 3 of the License, or
%%% any later version.
%%% 
%%% This program is distributed in the hope that it will be useful,
%%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%% GNU General Public License for more details.
if nargin<2 || isempty(thresh)
   thresh = 150; 
end
 
[M V] = openIMG(fn);
ss = size(M);

M2 = reshape(M, prod(ss(1:3)), size(M,4))';
avg = nanmean(M2);
sd = nanstd(M2);
snr = avg./sd;

AVG = reshape(avg', size(M,1), size(M,2), size(M,3));
SD = reshape(sd', size(M,1), size(M,2), size(M,3));
SNR = reshape(snr', size(M,1), size(M,2), size(M,3));

tmp = find(fn==filesep);
if isempty(tmp)
    pp = [pwd filesep]; 
    name = fn;
else
    pp = fn(1:tmp(end));
    name = fn(tmp(end)+1:end);
end

if ~exist([pp 'SNR_Images'],'dir')
    mkdir([pp 'SNR_Images']);
end

vv = V(1);
vv.fname = [pp 'SNR_Images/' name(1:end-4) '_MEANimage.nii'];
spm_write_vol(vv,AVG);
vv.fname = [pp 'SNR_Images/' name(1:end-4) '_SDimage.nii'];
spm_write_vol(vv,SD);
vv.fname = [pp 'SNR_Images/' name(1:end-4) '_SNRimage.nii'];
spm_write_vol(vv,SNR);

ind = find(avg>thresh);

S.Label = 'Global';
S.AVG = nanmean(avg(ind));
S.SD  = std(nanmean(M2(:,ind),2));
S.SNR = S.AVG./S.SD;


if nargin==3
   for ii = 1:length(ROIs);
       [mm hh] = openIMG(ROIs{ii});
       mm = resizeVol2(hh,V(1));
       ind2 = find(mm>0);
       ind3 = intersect(ind,ind2);
       
       lab = ROIs{ii};
       in1 = find(lab==filesep);
       in2 = find(lab=='.');
       lab = lab(in1(end)+1:in2(end)-1);
       
       S(ii+1).Label = lab;
       S(ii+1).AVG = nanmean(avg(ind3));
       S(ii+1).SD = nanstd(nanmean(M2(:,ind3),2));
       S(ii+1).SNR = S(ii+1).AVG./S(ii+1).SD;
   end
end

WriteDataToText(S,[pp 'SNR_Images/Global_SNR_Report_' name(1:end-4) '.txt'],'w','\t');
type([pp 'SNR_Images/Global_SNR_Report_' name(1:end-4) '.txt']);

snr = nanmean(M2).*nanstd(M2);
snr(find(snr<1))=0;
vv.fname = [pp 'SNR_Images/' name(1:end-4) '_aps_sig.nii'];
vv.dt = [64 0];
NV = zeros(ss(1:3))*NaN;
NV(:) = snr;
spm_write_vol(vv,NV);

