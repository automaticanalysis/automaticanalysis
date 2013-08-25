% Transform a 3D RGB (i.e. 4D data) to format displayable with montage
% Tibor Auer MRC CBU Cambridge 2012-2013

function mon = tr_RGBtoMontage(RGB)
for s = 1:size(RGB,3)
    mon(:,:,:,s) = squeeze(RGB(:,:,s,:));
end