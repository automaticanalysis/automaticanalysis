%% Automatic analysis - setup qsub
% To set up torque (qsub) and its relevant time and memory estimates, run
% this function when setting up study specific options in your aa_user
% script
% Also kills previous torque jobs that may have been running.
% matrixSize: either a 3x1 EPI matrix dimension vector or a single voxel count
% numImages: number of EPI images in experiment

function aap = aas_setup_qsub(aap, matrixSize, numImages)

% Default values for 2D EPI
if nargin < 3
    numImages = 2000;
end
if nargin < 2
    matrixSize = [64 64 32];
end

if length(matrixSize) ~= 1
    matrixSize =  matrixSize(1)*matrixSize(2)*matrixSize(3);    
end

% Add the path with functions to interact with torque (qsub) <-- changed
addpath(genpath(fullfile(aap.directory_conventions.fieldtripdir, 'qsub')))

% Remove previous jobs
qsublist('killall')

% Process in torque (qsub) rather than locally
aap.options.wheretoprocess = 'qsub'; 
aap.options.timelog = 0; % Already timelogged by torque (qsub)

% The timeBase and memoryBase estimates of the modules are based on typical
% EPI data of 64x64x32 matrix size, and on ~2000 EPI images (in 2 sessions)

% Assume that memory changes linearly with number of voxels in study
aap.options.qsub.memoryMult = matrixSize ./ (64 * 64 * 32) ...
    * numImages ./ 2000;

% Assume that time changes linearly with number of voxels in study
aap.options.qsub.timeMult = matrixSize ./ (64 * 64 * 32) ...
    * numImages ./ 2000;
