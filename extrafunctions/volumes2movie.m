function volumes2movie(V,outfile,fps)
% Converts a set of spm volumes (V struct) to a movie, saved as outfile and
% played back at fps.
% volumes2movie(V,outfile,fps)
xyz = spm_read_vols(V);
% Take 3 axial sections
steps = round(size(xyz,3)/4);
secs = [steps steps*2 steps*3];
% Flip around to appropriate orientation
for s = 1:3
    axials{s} = permute(squeeze(xyz(:,:,secs(s),:)),[2 1 3]);
    axials{s} = axials{s}(end:-1:1,:,:);
end
% And concatenate
axials = [axials{1} axials{2} axials{3}];
% And get rid of ridiculously big variable
clear xyz

% set video in standard display format to avoid weird distortions in most
% video players
axsize = size(axials);
goodsizes = [240 320; 480 640; 600 800; 768 1024];
ind = find(axsize(2)<goodsizes(:,2),1,'first');
targsize = goodsizes(ind,:);
if axsize(1)<targsize(1)
    % add 0s at bottom
    axials((axsize(1)+1):targsize(1),:,:) = 0;
end
if axsize(2)<targsize(2)
    axials(:,(axsize(2)+1):targsize(2),:) = 0;
end

% Rescale and make uint8
intmin = min(axials(:));
axials = axials - intmin;
intmax = max(axials(:));
axials = (axials/intmax);
axials = gray2ind(axials,255);

% Make movie frames
for f = 1:size(axials,3)
    mframes(f,:) = im2frame(cat(3,axials(:,:,f),axials(:,:,f),axials(:,:,f)));
end
% Write out avi
movie2avi(mframes,outfile,'FPS',fps,'COMPRESSION','none');
