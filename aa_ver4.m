function aa_ver4
% Automatic analysis - adds paths for aa commands to Matlab path list

aa_ver4_nocloud;

aapath = fileparts(mfilename('fullpath'));
fprintf('Please wait moment, adding Java objects\n');
run(fullfile(aapath,'cloudclient','aacloudclient.m'));
