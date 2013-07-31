function aa_ver4_nocloud
% Automatic analysis - adds paths for aa commands to Matlab path list

aapath = fileparts(mfilename('fullpath'));
addpath(genpath(aapath)); % recursively add AA subfolders

[pth nme ext]=fileparts(aapath);

fprintf('Welcome to aa version 4.1 on github July 2012\n',nme);
