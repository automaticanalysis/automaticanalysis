function mask_img(Mimg, fns, maskVal)
if nargin < 3 || isempty(maskVal)
    maskVal = NaN;
end
if ischar(Mimg)
    % Load mask
    M = spm_read_vols(spm_vol(Mimg));
end
if ischar(fns)
    fns = strvcat2cell(fns);
end    

for f = 1:length(fns)
    % Load image
    V = spm_vol(fns{f});
    Y = spm_read_vols(V);
    
    if isempty(Mimg)
        % Set things that are 0 to NaN
        Y(or(Y==0, ~isfinite(Y))) = maskVal;
    elseif all(size(M) == size(Y))
        % Mask image
        Y(~M) = maskVal;
    elseif any(size(M) ~= size(Y))
        aas_log([],true,'Mask and image to be masked are not of the same size!')
    end
    
    % Write image back...
    spm_write_vol(V,Y);
end