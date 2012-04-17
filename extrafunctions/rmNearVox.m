% This function erodes an image using another:
% The eroded image (mask2erode) is checked against another (erodingMask)
% so that any voxels in "mask2erode" at a certain radius from voxels in
% "erodingMask" are removed.
function mask2erode = rmNearVox(mask2erode, erodingMask, minDist)

% Checks if our minDist is an integer (necessary for smoothing function!)
if mod(minDist,1)==0
    % Then try the blunt approach
    erodingMask = smooth3(erodingMask, 'box', [minDist*2+1 minDist*2+1 minDist*2+1]);
    mask2erode(logical(erodingMask)) = 0;
else
    error('minDist must be an odd integer, as smoothing function will not accep otherwise')
end