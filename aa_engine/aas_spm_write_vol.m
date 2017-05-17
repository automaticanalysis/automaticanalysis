% function aas_spm_write_vol(Vout,Y,outfn)
% Like spm_write_vol, but handles 3D or 4D 
% For 4D, Vout can be EITHER structure same size as Y 
%              OR SINGLE structure 
%     In either case, timepoints dimension is set to 1 for each image plane
% Y can be 3D or 4D
% If optional third argument outfn is provided, the filename is set to this
% Rhodri Cusack BMI Western 2016-08-29

function aas_spm_write_vol(Vout,Y,outfn)

nimg=size(Y,4);

if length(Vout)==1
    Vout=repmat(Vout,[nimg 1]);
end;

for imgind=1:nimg
    if exist('outfn','var')
        Vout(imgind).fname=outfn;
    end;
    Vout(imgind).dim=Vout(imgind).dim(1:3);
    Vout(imgind).n=[imgind 1];
    spm_write_vol(Vout(imgind),Y(:,:,:,imgind));
end