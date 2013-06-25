function aa_ver4
% Automatic analysis - adds paths for aa commands to Matlab path list

aapath = fileparts(mfilename('fullpath'));
addpath(genpath(aapath)); % recursively add AA subfolders

[pth nme ext]=fileparts(aapath);

fprintf('Welcome to aa version 4 on github May 2012\n',nme);
fprintf('Please wait moment, adding Java objects\n',nme);
run(fullfile(aapath,'cloudclient','aacloudclient.m'));
