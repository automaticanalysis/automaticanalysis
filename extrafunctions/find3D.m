function loc = find3D(dat)
if ndims(dat) == 3
    for i = 1:size(dat,3), if ~isempty(find(dat(:,:,i), 1)), break; end; end
    loc(3) = i;
    dat = dat(:,:,i);
end

for i = 1:size(dat,2), if ~isempty(find(dat(:,i), 1)), break; end; end
loc(2) = i;
dat = dat(:,i);

loc(1) = find(dat,1);

