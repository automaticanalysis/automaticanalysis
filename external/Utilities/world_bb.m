function [bb tc] = world_bb(V)
%%%  world-bb -- get bounding box in world (mm) coordinates
%%%  I scooped this from John Ashburner's Website
%%%
%%% Inputs:
%%% V = a volume header (e.g. V = spm_vol(fn);)
%%%
%%% Outputs:
%%%
%%% bb = the bounding box
%%% tc = the locations of all 6 corners in world space
%%%
%%%
%%% From http://www.cs.ucl.ac.uk/staff/g.ridgway/vbm/world_bb.m
%%% Based on John Ashburner's reorient.m
%%% http://www.sph.umich.edu/~nichols/JohnsGems.html#Gem7
%%% http://www.sph.umich.edu/~nichols/JohnsGems5.html#Gem2
%%% Adapted by Ged Ridgway -- email bugs to drc.spm@gmail.com
%%%
%%% Reeditied by Aaron Schultz - aschultz@martinos.org

d = V.dim(1:3);
% corners in voxel-space
c = [ 1    1    1    1
    1    1    d(3) 1
    1    d(2) 1    1
    1    d(2) d(3) 1
    d(1) 1    1    1
    d(1) 1    d(3) 1
    d(1) d(2) 1    1
    d(1) d(2) d(3) 1 ]';
% corners in world-space
tc = V.mat(1:4,1:4)*c;
tc = tc(1:3,:);

% bounding box (world) min and max
mn = min(tc,[],2)';
mx = max(tc,[],2)';
bb = [mn; mx];
