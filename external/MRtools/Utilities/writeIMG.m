function writeIMG(V,M,nfn,bitdepth)
%%%
%%% This a simple wrapper for to write spm and freesurfer image data to file.
%%% You will need an spm and freesurfer installation respectively for this
%%% function to work.
%%%
%%% Inputs:
%%% V = a volume header
%%% M = a 2D, 3D, or 4D matrix of image data
%%% nfn = output filename
%%%
%%%
%%% Written by Aaron Schultz - aschultz@martinos.org
%%% Copyright (C) 2011,  Aaron P. Schultz
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

if exist(nfn);
    warning off;
    delete(nfn);
    warning on;
end

V = V(1);

if isfield(V,'vox2ras0')
    V.descrip = [];
    V.vol = reshape(M', V.ss);
    
    if nargin == 4
        MRIwrite(V,nfn,bitdepth);
    else
        MRIwrite(V,nfn,'float');
    end
else
    for ii = 1:size(M,4);
        if numel(V)>1
            V(ii).fname = nfn;
            V(ii).n = [ii 1];
            
            if nargin==4
                V(ii).dt = [bitdepth 0];
            end
            spm_write_vol(V(ii),M(:,:,:,ii));
        else
            V.fname = nfn;
            V.n = [ii 1];
            
            if nargin==4
                V.dt = [bitdepth 0];
            end
            spm_write_vol(V,M(:,:,:,ii));
        end
    end
end




