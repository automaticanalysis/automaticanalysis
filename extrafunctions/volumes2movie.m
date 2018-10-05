function volumes2movie(V,outfile,fps)
%
% Convert a set of spm volumes (V struct) into an avi movie, saved as 
% outfile and played back at fps. Movie is 3x2 axial sections that span
% the entire z volume of the brain, each showing voxels averaged 
% across the substack (ergo "glassbrain" movie)
%
% NB: Orientation is along positive z (i.e. bottom-to-top). Think: standing 
% at subject's feet looking into the bore of the scanner, so in each section
% face is up but right and left are reversed. 
%
% CHANGE HISTORY
%
% 1/2018 [MSJ]	new, based on volume2movie
%

xyz = spm_read_vols(V);

% generate 6 axial sections that span stack

steps = round(size(xyz,3)/7);
secs = [ 1 steps steps*2 steps*3 steps*4 steps*5 steps*6 ];

% Flip around to appropriate orientation, averaging across each substack

for s = 1:6
    axials{s} = permute(squeeze(mean(xyz(:,:,secs(s):secs(s+1),:),3)),[2 1 3]);
    axials{s} = axials{s}(end:-1:1,:,:);
end

% concatenate to 3x2

axials = [axials{1} axials{2} axials{3} ; axials{4} axials{5} axials{6}];

%  get rid of mem hog variable

clear xyz

% pick a standard display format to avoid distortions 

axsize = size(axials);
goodsizes = [240 320; 480 640; 600 800; 768 1024];
ind = find(axsize(2)<goodsizes(:,2),1,'first');
targsize = goodsizes(ind,:);

% pad to fit

if axsize(1)<targsize(1)
    axials((axsize(1)+1):targsize(1),:,:) = 0;
end

if axsize(2)<targsize(2)
    axials(:,(axsize(2)+1):targsize(2),:) = 0;
end

% rescale and make uint8

intmin = min(axials(:));
axials = axials - intmin;
intmax = max(axials(:));
axials = (axials/intmax);
axials = gray2ind(axials,255);

% quick-n-dirty contrast cleanup
% however, this requires image processing toolbox

if (exist('imadjust'))
	for index = 1:size(axials,3)
		axials(:,:,index) = imadjust(axials(:,:,index));
	end
end

% make movie frames

for f = 1:size(axials,3)
    mframes(f,:) = im2frame(cat(3,axials(:,:,f),axials(:,:,f),axials(:,:,f)));
end

% write out avi

% movie2avi(mframes,outfile,'FPS',fps,'COMPRESSION','none');
%
% movie2avi is no longer available. Use writeVideo instead:

v = VideoWriter(outfile, 'Uncompressed AVI');
v.FrameRate = fps;
open(v)
writeVideo(v, mframes);
close(v);
