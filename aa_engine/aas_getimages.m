% Automatic Analysis function to retrieve image lists
% function [imagefns]=aas_getimages(aap,i,j,prefixes,minimg,maximg)
%  select images from subject i session j
%  only images beginning with 'prefixes' are returned
%  either provide range using minimg, maximg (0,inf for all) or provide
%   single parameter to get one image
% Rhodri Cusack MRC CBU Cambridge, Nov 2005

function [imagefns]=aas_getimages(aap,i,j,prefixes,minimg,maximg)

if (~exist('minimg','var'))
    minimg=0;
    maximg=inf;
else
    if (~exist('maximg','var')), maximg=minimg; end;
end;

sesspath=fullfile(aas_getsubjpath(aap,i),aap.acq_details.sessions(j).name);

if (aap.directory_conventions.selectechonumbers)
    allimages=dir(fullfile(sesspath,sprintf('%s%s*-%02d.nii',prefixes,aap.directory_conventions.rawdataafterconversionprefix,aap.acq_details.sessions_echonumbers(j))));
else
    allimages=dir(fullfile(sesspath,sprintf('%s%s*.nii',prefixes,aap.directory_conventions.rawdataafterconversionprefix)));
end;

imagefns=[];
pth=sesspath;
fns=[];
for i=1:length(allimages)
    [t r]=strtok(allimages(i).name,'-'); % start of series num
    [t r]=strtok(r,'-'); % start of imgnum
    t=strtok(r,'-');
    imgnum=str2num(t);
    if (imgnum>=minimg & imgnum<=maximg)
        imagefns=strvcat(imagefns,fullfile(sesspath,allimages(i).name));
        fns=strvcat(fns,allimages(i).name);    
    end;
end;
