%% This function checks the SPM model of your choice at the 1st level
% Inputs:
%   pth = path to your model (defaults to current path)
%   sess = which sessions to examine (all by default
%   maskVol = which mask to use to mask in your voxels of interest
%       (defaults to using all voxels that are finite and not 0);
function h = firstlevelmodelStats(pth, sess, maskVol)
if nargin < 1
    pth = pwd;
end
if nargin < 2
    sess = [];
end
if nargin < 3
    maskVol = [];
end

h = [];

% Get betas, SPM, maskVol
D = dir(fullfile(pth, 'beta_*'));
load(fullfile(pth, 'SPM.mat'));
SPMcols = SPM.xX.iC;
if ~isempty(sess)
    SPMsess = [];
    for s = sess
        SPMsess = [SPMsess, SPM.Sess(s).col];
    end
    SPMcols = SPMcols(ismember(SPMcols, SPMsess));
end
SPMnames = SPM.xX.name;
if ~isempty(maskVol)
    if ischar(maskVol)
        maskVol = spm_read_vols(spm_vol(maskVol));
    end
    maskVol = logical(maskVol);
end

%% INTEREST AND NUISANCE COLUMNS...
SPMinterest = SPM.xX.name(SPM.xX.iC);
SPMnuisance = SPM.xX.name(SPM.xX.iG);

errorCols = [];

% Detect common errors...
fprintf('Columns of interest  \n');
for n = 1:length(SPMinterest)
    fprintf([SPMinterest{n}  '\n'])
    tmp = SPMinterest{n}(7:end);
    if ~isempty(strfind(tmp, 'Spike')) || ... % Spikes
            ~isempty(strfind('xyzrpj', tmp)) || ... % Movement pars
            strcmp('GM', tmp) || strcmp('WM', tmp) ||  strcmp('CSF', tmp) ||  strcmp('OOH', tmp) % CompRegs
        errorCols = [errorCols SPM.xX.iC(n)];
    end
end
fprintf('\n');
fprintf('Columns of nuisance \n');
for n = 1:length(SPMnuisance)
    fprintf(['\t' SPMnuisance{n} '\n'])
end

if ~isempty(errorCols)
   fprintf('\n')
   fprintf('Problems with columns... \n')
   for n = 1:length(errorCols)
       fprintf(['\t' SPMnames{errorCols(n)}  '\n'])
   end
end
if ~isempty(errorCols)
    error('Error in your model!')
end

%% CORRELATION BETWEEN REGRESSORS...
SPMmodel = SPM.xX.X(:, SPM.xX.iC);
[sharedVar, h.regs] = corrTCs(SPMmodel, SPMnames(SPM.xX.iC), 1, 0);

%% CORRELATIONS BETWEEN BETAS...
fprintf('Loading beta images into data structure\n')
for d = 1:length(SPMcols)
    V = spm_vol(fullfile(pth, D(SPMcols(d)).name));
    Y = spm_read_vols(V);
    % maskVol things we don't want...
    if ~isempty(maskVol);
        Y = Y(maskVol);
    else
        Y = Y(isfinite(Y) & Y~=0);  
    end
    
    if d == 1
        data = nan(length(Y), length(SPMcols));
    end    
    data(:,d) = Y;
end

fprintf('Correlating voxel-wise betas across regressors\n')
[sharedVar, h.betas] = corrTCs(data, SPMnames(SPMcols), 1, 0);
