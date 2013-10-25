function fake_SPM(spmT_images, con_images, df, analysis_fn)
%% fake_SPM(spmT_images, con_images, df, analysis_fn)
% input:
% 1) spmT_image(s): cell array of image names containng t-values
% 2) con_image(s): cell array of image names containing effect/intensity
%      information (e.g. contrast value, correlation index, etc.)
% 3) df: single value (used for all images) or vector (for respective images)
%       containing degrees of freedom for statistical contrast
% 4) analysis_fn: the analysis folder name, where output will be saved
%
% NOTE: all the images *MUST* be in the same image space
%
% output: Will save spmT_%04d.img, con_%04d.img and SPM.mat

%% Check inputs
if nargin < 4
    analysis_fn = 'fake_SPM';
end
if nargin < 3
    error('You need to provide the spmT_images, the con_images and the df!')
end
% Make image inputs into cells if they are not (@@@ IMPROVE LATER @@@)
if ~iscell(spmT_images)
    spmT_images = strvcat2cell(spmT_images);
end
if ~iscell(con_images)
    con_images = strvcat2cell(con_images);
end
% Get df out of cell array
if iscell(df)
    %df = cell2mat(df);
    error('Degrees of freedom must be a single value for all images')
end
% Make sure equal number of images...
if any(size(spmT_images) ~= size(con_images))
    error('The number of spmT_images and con_images is different!')
end
%{
% Make degrees of freedom length of number of images
if any(size(spmT_images) ~= size(df))
    if all(size(df) == 1)
        df = repmat(df, size(spmT_images));
    else
        error('The number of spmT_images and df is different!')
    end
end
%}

%% Copy all our images into the analysis folder
if ~exist(analysis_fn, 'dir')
    mkdir(analysis_fn);
end
for c = 1:length(spmT_images)
    V = spm_vol(con_images{c});
    Y = spm_read_vols(V);
    V.fname = fullfile(analysis_fn, sprintf('con_%04d.nii', c));
    spm_write_vol(V, Y);
    
    V = spm_vol(spmT_images{c});
    Y = spm_read_vols(V);
    V.fname = fullfile(analysis_fn, sprintf('spmT_%04d.nii', c));
    spm_write_vol(V, Y);
end

%% Create empty SPM structure!
SPM = [];
SPM.xX.xKXs.X = 1;
SPM.xX.xKXs.tol = 1;
SPM.xX.xKXs.ds = 1;
SPM.xX.xKXs.u = 1;
SPM.xX.xKXs.v = 1;
SPM.xX.xKXs.rk = 1;
SPM.xX.xKXs.oP = [];
SPM.xX.xKXs.oPp = [];
SPM.xX.xKXs.ups = [];
SPM.xX.xKXs.sus = [];

SPM.xX.name = {'fake_SPM'};
SPM.xX.erdf = df;
SPM.xX.X = 0;
SPM.xX.nKX = 0;

SPM.xY.VY = [];
SPM.xY.VY.fname = {'fake_SPM'};

SPM.VResMS = [];
SPM.Vbeta = [];


% Clear SPM.xCon
SPM.xCon = [];

% Set correct path
SPM.swd = analysis_fn;

% Set world coordinates for visualisation...
% ...which should already be found in the images...
SPM.xVol.M = V.mat;
SPM.xVol.iM = inv(SPM.xVol.M);

% Size of the volume
SPM.xVol.DIM = V.dim';

% Smoothness of the volume...
% ...Get the number of mm per voxel...
mmVox = vox2mm(V);

% ...then set the FWHM
% @@@ FIND WAY TO ESTIMATE FWHMmm... @@@
FWHMmm = min(mmVox./2);

SPM.xVol.FWHM = [FWHMmm FWHMmm FWHMmm];
SPM.xVol.FWHM = SPM.xVol.FWHM ./ mmVox;

% Spm_resels_vol function
% NOTE: This is probably not valid for FWE still, since the
% searchlight procedure means each voxels is already "smoothed" to
% some extent...
SPM.xVol.R = spm_resels_vol( ...
    spm_vol(fullfile(analysis_fn, 'con_0001.nii')), ...
    SPM.xVol.FWHM)';

% Included voxels
[X Y Z] = ind2sub(SPM.xVol.DIM', find(and(Y~=0, isfinite(Y))));
SPM.xVol.XYZ = [X';Y';Z'];

% Length of voxels in analysis
SPM.xVol.S = length(X);

% Filehandle of resels per voxel image (i.e. none!)
SPM.xVol.VRpv = [];

for c = 1:length(spmT_images)
    % SPM.xCon (.name)
    SPM.xCon(c).name = sprintf('con%d', c);
    SPM.xCon(c).STAT = 'T';
    SPM.xCon(c).c = 1; % SPM.xCon(c).c = ones(df+1,1); %ones(df(c)+1,1);
    SPM.xCon(c).eidf = df; %df(c);
    SPM.xCon(c).Vcon = spm_vol(fullfile(analysis_fn, sprintf('con_%04d.nii', c)));
    SPM.xCon(c).Vspm = spm_vol(fullfile(analysis_fn, sprintf('spmT_%04d.nii', c)));
end

% Save SPM
save(fullfile(analysis_fn, 'SPM.mat'), 'SPM');

end