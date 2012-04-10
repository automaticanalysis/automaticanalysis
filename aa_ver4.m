% Automatic analysis - adds paths for aa commands to Matlab path list
% If you are installing aa you must change the first line of this script

aapath='/home/rhodricusack/software/automaticanalysis/';
addpath(fullfile(aapath,'aa_engine'));
addpath(fullfile(aapath,'aa_modules'));
addpath(fullfile(aapath,'aa_recipes_and_parametersets'));
addpath(fullfile(aapath,'aa_config'));
addpath(genpath(fullfile(aapath,'aa_toolbox')));
addpath(fullfile(aapath,'cbusoftware'));
addpath(fullfile(aapath,'examples'));
addpath(fullfile(aapath,'extrafunctions'));
[pth nme ext]=fileparts(aapath);
fprintf('Welcome to aa version 4 on github March 2012\n',nme);
fprintf('Please wait moment, adding Java objects\n',nme);
run(fullfile(aapath,'cloudclient','aacloudclient.m'));
