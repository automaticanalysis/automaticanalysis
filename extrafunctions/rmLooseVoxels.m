% Function rmLooseVox
% Inputs:
% imgM - image Matrix, or file name
% connThres - connectome threshold (when does some bit is classified as a
% loose bit)

function imgM = rmLooseVoxels(imgM, connThresh)

if varagin < 2
    connThresh = 6;
end

if ischar(imgM)
    [Mpth, Mfn, Mext] = fileparts(imgM);
    V = spm_vol(imgM);
    Y = spm_read_vols(V);
elseif isnumeric(imgM);
    Y = imgM;
else
    error('Not a valid input. Input an image filename or image matrix.')
end

CC = bwconncomp(Y, connThresh);
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,idx] = max(numPixels);
numPixels(idx) = []; % For display purposes...
sumPixels = sum(numPixels);

if length(CC.PixelIdxList) == 1
else
    for b = 1:CC.NumObjects
        if b~=idx
            Y(CC.PixelIdxList{b}) = 0;
        end
    end
end
fprintf('The brain mask consists of %d clusters, of approx. %0.2f voxels (main cluster = %d)\n', length(CC.PixelIdxList), mean(numPixels), biggest)

if ischar(imgm)
    V.fname = fullfile(mpth, ['l' mfn mext]);
    spm_write_vol(v,y);
    spm_check_registration(Mimg)
    spm_orthviews('addcolouredimage',1,V.fname, [1 0 0])
end

end
