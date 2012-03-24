% epi can be filename or volume
% function [oftemp]=aas_overlapwithtemplate(aap,epifn,typeofthresh)
% typeofthresh=1 - use 10% of distance between 3 percentile & 97 percentile
% typeofthresh=2 - everything but 0 or nans


function [oftemp]=aas_overlapwithtemplate(aap,epifn,typeofthresh)

try
    Vtemplate
catch
    Vtemplate=spm_vol(aap.directory_conventions.T1template);
end;

[Ytemplate XYZ]=spm_read_vols(Vtemplate);

XYZ=[XYZ;ones(1,size(XYZ,2))];

Vepi=spm_vol(epifn);
XYZ=Vepi.mat\XYZ;

Yepi=spm_sample_vol(Vepi,XYZ(1,:),XYZ(2,:),XYZ(3,:),0);
Yepi=reshape(Yepi,size(Ytemplate));

% thresholda (proportion of distance between 3 percentile & 97 percentile -
% used in image processing as a measure robust to outliers)

threshprop=0.10;

% work out threshold for template
ys=sort(Ytemplate(~isnan(Ytemplate)));
bright3=ys(round(length(ys)*0.3));
bright97=ys(round(length(ys)*0.97));
thresh=bright3*(1-threshprop)+bright97*threshprop;
Ytemplate=Ytemplate>thresh;

if (typeofthresh==1)
    % work out threshold for epi
    ys=sort(Yepi(~isnan(Yepi)));
    bright3=ys(round(length(ys)*0.3));
    bright97=ys(round(length(ys)*0.97));
    threshepi=bright3*(1-threshprop)+bright97*threshprop;
    Yepi=Yepi>threshepi;
elseif (typeofthresh==2)
    Yepi(isnan(Yepi))=0;
    Yepi=~(Yepi==0);
end;

% figure(10);
% subplot(211);
% imagesc(Yepi(:,:,20));
% subplot(212);
% imagesc(Ytemplate(:,:,20));

Yepi=Yepi(:);
Ytemplate=Ytemplate(:);


overlap=sum(Yepi.*Ytemplate);
oftemp=overlap/sum(Ytemplate);





