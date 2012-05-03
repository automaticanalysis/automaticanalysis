% Automatic analysis - setup qsub

function aap = aas_setup_qsub(aap, matrixSize, numImages)

% Add qsub path
addpath(aap.directory_conventions.fieldtripdir)

% Remove previous workers
qsublist('killall')

aap.options.wheretoprocess = 'qsub';
aap.options.timelog = 0; % Already timelogged by qsub

% The timeBase and memoryBase estimates of the modules are based on typical
% EPI data of 64x64x32 matrix size, and on ~2000 EPI images (in 2 sessions)

% Assume that memory changes linearly with number of voxels in study
aap.options.qsub.memoryMult = matrixSize(1)*matrixSize(2)*matrixSize(3) ./ (64 * 64 * 32) ...
    * numImages ./ 2000;

% Assume that time changes linearly with number of voxels in study
aap.options.qsub.timeMult = matrixSize(1)*matrixSize(2)*matrixSize(3) ./ (64 * 64 * 32) ...
    * numImages ./ 2000;