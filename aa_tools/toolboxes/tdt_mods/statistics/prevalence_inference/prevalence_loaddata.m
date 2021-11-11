% function [a, mask, vol] = prevalence_loaddata(fnames_in, decoding_measure, mask)
%
% Helper function to load data for the prevalence analysis and prepare some
% variables. Can handle SPM images (.nii, .img) and TDT result .mat files
% (searchlight, ROI, wholebrain).
%
% IN
%   fnames_in: nsbj x n_images cellstr, where fnames_in(:, 1)
%       contains the original unpermuted image.
%       Supported image types:
%          * .img, .nii: image files that can be read with SPM
%          *       .mat: resultfile from TDT (Todo CHARACTERISTICS)
%       Images need to be spatially coregistered and resliced, so that the
%       same voxel in each image corresponds to the same location in all
%       other images.
% OPTIONAL
%   decoding_measure (TDT only): String with decoding_measure to select
%       from result (e.g. 'accuracy_minus_chance'). Only necessary if the
%       result file contains more than one decoding_measure.
%   mask: Only for SPM at the moment: Logical vector with inmask = true
%       voxels, if the mask should not be determined from the data. Should
%       have the dimensions of the original images.
% OUT
%   a: ndim x N x P1 data matrix, where
%      ndim: number of inmask voxels (or in general number of dimensions,
%           e.g. ROIs)
%         N: number of subjects (or in general units of observations)
%        P1: number of permutations per subject (unit)
%       a(:, :, 1) should contain all original unpermuted result images of
%          all subjects that will be used for testing.
%       a(:, k, i) contains all inmask data from fnames_in{k, i};
%   mask: 1d/2d/3d logical matrix that contains the mask used to truncate
%       the data. size(mask) is the size of the original input images.
%       The mask is typically retrieved from the data.  
%       All voxels that are in all input images always nan or 0 are set to
%       be out of the mask (value set to false).
%       The value of nvox is always larger than ndim above, because 
%       ndim = sum(mask(:)==true).
%   vol: Will contain a .mat field with the 4x4 rotation and translation
%       matrix. If this cannot be extracted from the input image, the field
%       is left empty. If SPM is used, vol will contain the full volume
%       information of the first image. If TDT is used, .datainfo from the
%       result will be available.
%
% Kai, 2016/07/25

function [a, mask, vol] = prevalence_loaddata(fnames_in, decoding_measure, mask)

%% Determine what type of data we have
if ~iscellstr(fnames_in)
    error('prevalence_loaddata:invalid_fnames_in', 'Please provide a cellstr with files to load')
end

% get extension of first input file to determine file type
[fd, fn, ext] = fileparts(fnames_in{1});
if  strcmp(ext, '.mat')
    inputformat = 'mat';
    disp = ['Loading data from ' ext ' files'];
elseif strcmp(ext, '.img') || strcmp(ext, '.nii') || strcmp(ext, '.hdr') || strcmp(ext, '.gz')
    inputformat = 'SPM';
    disp = ['Loading data from ' ext ' files using SPM'];
else
    error('prevalence_loaddata:unkown_dataformat', 'Unkown extension "%s", please implement detection and loading', ext);
end
fprintf('Datatype has been determined as: %s\n', inputformat);

%% Load all data
[N, P1] = size(fnames_in);
fprintf('loading data\n')
% load accuracy images
a = cell(N, P1);

%% Check if mask is provided, and handle it
if exist('mask', 'var')
    if ischar(mask)
        error('Loading an external mask from file is not implemented at the moment, consider to implement it. If you do this, please check the orientation below')
    elseif islogical(mask)
        warning('Taking external mask only works for SPM at the moment, but even there we cannot check the orientation (will abort at the end if the number of voxels is not equal between mask and images).')
    else
        error('The external mask only works for SPM at the moment, and needs to be provided as logical vector');
    end
end

%% For all subjects
for k = 1 : N
    fprintf('  subject #%d: ', k)
    
    %% Load all images
    for i = 1 : P1
        
        if isempty(fnames_in{k, i})
            error('prevalence_loaddata:empty_filename', 'No filename found in fnames_in for sbj %i, image %i. All subjects currently need the same amount of images. Please check or change code to deal with it.', k, i);
        end
        
        %% Read SPM
        if strcmp(inputformat, 'SPM')
            curr_fname = fnames_in{k, i};
            gz = strcmp(ext, '.gz'); % gzipped image?
            if gz
                curr_fname = gunzip(curr_fname, tempdir);   % gunzip to temporary directory
                curr_fname = curr_fname{1};
            end
            vol = spm_vol(curr_fname);
            Y = spm_read_vols(vol);
            if ndims(Y) == 4 % 4D input
                a(k, i:i+size(Y,4)-1) = arrayfun(@(v) reshape(Y(:,:,:,v),[],1), 1:size(Y,4), 'UniformOutput', false);
                vol = vol(1);
            else
                a{k, i} = Y(:);
            end
            if gz
                delete(curr_fname); % delete temporary file again
            end
            
            %% Read TDT
        elseif strcmp(inputformat, 'mat')
            %% Load TDT mat file
            currmat = load(fnames_in{k, i});
            
            %% Determine at first call which decoding measure to use if not provided as argument
            if ~exist('decoding_measure', 'var') || isempty(decoding_measure)
                % check if only 1 measure exists and if so use it
                fnames = fieldnames(currmat.results);
                decoding_measure_ind = structfun(@(x) isfield(x, 'output'), currmat.results); % find all fields that contain .output
                if sum(decoding_measure_ind) == 0
                    error('%s: No decoding measure has been provided and could not detect any decoding measure, please check.', fnames_in{k, i})
                elseif sum(decoding_measure_ind) >= 2
                    found_decoding_measures = fnames(decoding_measure_ind) %#ok<NOPRT,NASGU>
                    error('%s: Found multiple decoding measures (above), please provide one of those as argument to the function', fnames_in{k, i})
                else
                    decoding_measure = fnames{decoding_measure_ind};
                    dispv(1, '%s: Detected decoding measure: %s, using this for prevalence analysis', fnames_in{k, i}, decoding_measure);
                end
                clear fnames decoding_measure_ind
            end
            
            %% Load searchlight
            if ~exist('first_analysis', 'var')
                first_analysis = currmat.results.analysis;
            end
            if ~strcmp(first_analysis, currmat.results.analysis)
                error('The first TDT analysis was of type %s, but the current analysis is of type %s. All analyses must have the same analysis. Aborting', first_analysis, currmat.results.analysis);
            end
            
            a{k, i} = currmat.results.(decoding_measure).output;
            % check output format & size
            if ~isnumeric(a{k, i})
                error('%s: Input measure %s is not numeric, aborting',  fnames_in{k, i}, decoding_measure)
            elseif (size(a{k, i}, 1) ~= 1 && size(a{k, i}, 2) ~= 1) || length(size(a{k, i})) > 2
                size(a{k, i})
                error('%s: Input measure %s is not a vector but of size above, aborting', fnames_in{k, i}, decoding_measure)
            elseif (k>1 || i>1) && length(a{k, i})~= length(a{1,1}) % check same number of entries for all input data by checking to the first
                error('%s: Data has a different length (%i) than the data of image 1 (%s: %i), aborting', fnames_in{k, i}, length(a{k, i}), fnames_in{1, 1}, length(a{1,1}))
            end
            
            %% get transformation matrix and datainfo from loaded file
            if ~exist('vol', 'var')
                try
                    vol.mat = currmat.results.datainfo.mat;
                catch
                    warning('Could not get transformation matrix from loaded mat file. This might be the case if the mat files have been created with an old version of TDT. In this case, you need to provide the orientation at the end of the function')
                    vol.mat = [];
                end
                try
                    vol.datainfo = currmat.results.datainfo;
                catch
                    warning('Could not get datainfo from loaded mat file. This might be the case if the mat files have been created with an old version of TDT. In this case, you need to provide the orientation at the end of the function')
                    vol.datainfo = [];
                end
                
                % Copy other fields that should be copied from results
                % this will be checked for consistency between loaded files
                % and returned in vol
                copy_fields_to_vol = {'mask_index', 'mask_index_each', 'roi_names'};
                for f_ind = 1:length(copy_fields_to_vol)
                    curr_field = copy_fields_to_vol{f_ind};
                    if isfield(currmat.results, curr_field)
                        vol.(curr_field) = currmat.results.(curr_field);
                    end
                end
            end
            
            %% unkown input type
        else
            error('prevalence_loaddata:unkown_inputformat_during_loading', 'Unkown extension "%s" during loading, please implement loading', inputformat)
        end
        
        % check that orientation matrices are the same
        if ~exist('mat', 'var')
            mat = vol.mat; % save orientation matrix from first image to compare to others
        end
        if ~isempty(vol.mat) || ~isempty(mat)
            mat_diff = abs(vol.mat(:)-mat(:));
            tolerance = 32*eps(max(vol.mat(:),mat(:)));
            if any(mat_diff > tolerance) % like isequal, but allows for rounding errors
                error('prevalence_loaddata:different_rotation_translation_mat', 'Rotation & translation matrix of image in file \n %s \n is different from rotation & translation matrix of the first file.\n The .mat entry defines rotation & translation of the image.\n That both differ means that at least one of both has been rotated.\n Please use reslicing (e.g. from SPM) to have all images in the same position.', fnames_in{k, i})
            end
        end
        
        % check other fields for consistency
        if exist('copy_fields_to_vol', 'var')
            for f_ind = 1:length(copy_fields_to_vol)
                if isfield(currmat.results, curr_field)
                    if ~isequal(vol.(curr_field), currmat.results.(curr_field))
                        disp(vol.(curr_field))
                        disp(currmat.results.(curr_field))
                        error(['prevalence_loaddatat:inconsistentField_' curr_field], ...
                            'results.%s differ between first loaded file (%s) and current file (%s), please check.', curr_field, fnames_in{1}, fnames_in{k, i});
                    end
                end
            end
        end
        
        fprintf('.')
    end
    fprintf('\n')
end
% a is now a cell array of size N x P1, where each cell contains voxel
% values in one column vector
a = cell2mat(reshape(a, [1, N, size(a,2)]));
% a is now a matrix of size
%  (number of voxels) x N (number of subjects) x P1 (number of permutations per subject)

%% get mask (SPM)
if strcmp(inputformat, 'SPM')
    if exist('mask', 'var')
        if ~islogical(mask)
            error('mask should be a logical vector here if provided, please check why that''s not the case');
        end
        disp('using externally provided mask to mask files, assuming it is a vector for the moment (if you implement loading the mask, put it at front and then check the orientation matrix)');
    else
        % determine mask from data; out-of-mask voxels may be NaN or 0
        mask = all(all(~isnan(a), 2), 3) & ~any(all(a == 0, 3), 2);
        mask = reshape(mask, size(Y,1:3));
    end
    
    % truncate data to in-mask voxels
    a = a(mask, :, :);
    
    %% get mask (TDT)
elseif strcmp(inputformat, 'mat')
    
    % we already have only the inmaskvoxels, so our job here is to create a
    % mask image that can later be used to generate the outputimages, and
    % not to mask all data like for SPM
    
    if exist('mask', 'var')
        error('External masking for TDT not implement yet, consider to simply mask images afterwards')
    end
    
    if strcmpi(currmat.results.analysis, 'ROI')
        % save mask_index_each instead of mask
        for m_ind = 1:length(currmat.results.mask_index_each)
            mask{m_ind} = false(currmat.results.datainfo.dim);
            mask{m_ind}(currmat.results.mask_index_each{m_ind}) = true;
        end
    elseif strcmpi(currmat.results.analysis, 'searchlight') || strcmpi(currmat.results.analysis, 'wholebrain')
        mask = false(currmat.results.datainfo.dim);
        mask(currmat.results.mask_index) = true;
    else
        error(['Unkown TDT results.analysis=' currmat.results.analysis]);
    end
    
    %% get mask (through error for unkown type)
else
    error('Unkown type %s to get mask, please implement')
end