function fo = img_flipud(fi)
for i = 1:size(fi, 3)
    fo(:,:,i) = flipud(fi(:,:,i));
end
