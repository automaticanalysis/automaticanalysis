% Automatic analysis - adds paths for aa commands to Matlab path list
% If you are installing aa you must change the first line of this script

aapath='/home/common/matlab/spm_batch/aa/release-4.0-beta/';
addpath(fullfile(aapath,'aa_engine'));
addpath(fullfile(aapath,'aa_modules'));
addpath(fullfile(aapath,'aa_recipes_and_parametersets'));
addpath(fullfile(aapath,'aa_config'));
addpath(fullfile(aapath,'cbusoftware'));
addpath(fullfile(aapath,'examples'));
[pth nme ext]=fileparts(aapath);
fprintf('Welcome to aa version 4-devel August 2010\n',nme);
fprintf('Please wait moment, adding Java objects\n',nme);
run(fullfile(aapath,'cloudclient','aacloudclient.m'));
