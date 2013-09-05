% Sets the datatype of an image...
function imgSetDT(Mimg, DT)

% Load image header
V = spm_vol(Mimg);
if V.dt(1) == DT
    fprintf('Image already has required datatype\n')
else
    % Get image matrix
    Y = spm_read_vols(V);
    
    % Set datatype...
    V.dt(1) = DT;
    
    % Write image back...
    spm_write_vol(V,Y);
end
