function aa_ver4_nocloud
% Automatic analysis - adds paths for aa commands to Matlab path list

aapath = fileparts(mfilename('fullpath'));
addpath(genpath(aapath)); % recursively add AA subfolders

% remove GitHub-related path
rmpath(genpath(fullfile(aapath,'.git')));

[pth nme ext]=fileparts(aapath);
fprintf('Welcome to aa version 4.2.0 Jan 2015\n');
fprintf(' If you publish work that has used aa, please cite our manuscript:\n');
fprintf(' Cusack R, Vicente-Grabovetsky A, Mitchell DJ, Wild CJ, Auer T, Linke AC, Peelle JE (2015) Automatic analysis (aa): Efficient neuroimaging workflows and parallel processing using Matlab and XML. Frontiers in Neuroinformatics 8:90.\n');
fprintf(' http://journal.frontiersin.org/Journal/10.3389/fninf.2014.00090/abstract\n');
