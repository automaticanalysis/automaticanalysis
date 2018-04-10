function [M h ind] = FastRead(fn,msk,ignore)
%%% Alternative tool for reading in NIFTI data.  Seems to perform a little
%%% better than openIMG.m
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

if isstruct(fn)
    h = fn;
else
    h = spm_vol(char(fn));
end


if ~any(h(1).dim==1)
    spm_check_orientations(h);
end

% if mean(mean(std(reshape([h.mat],4,4,numel(h)),[],3)))~=0
%     error('Images are not all the same orientation');
% end

I.v = h(1);
FullIndex = 1:prod(I.v.dim);

if nargin > 1
    if ischar(msk)
        mh = spm_vol(msk);
        mask = resizeVol(mh,h(1));
        FullIndex = find(mask==1);
    else
        FullIndex = msk;
    end
end

[x y z] = ind2sub(I.v.dim,FullIndex);
M = zeros(numel(h),numel(x));
% persisText();
for ii = 1:numel(h);
%     persisText([num2str(ii) ' of ' num2str(numel(h))]);
    M(ii,:) = spm_sample_vol(h(ii),x,y,z,0);
end
% persisText();


% toc

% tic
% m = openIMG(fn);
% m = reshapeWholeBrain(size(m),m);
% toc
