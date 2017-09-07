%% This function compares the output of normalised spmT distributions
%% across subjects
% Inputs:
%   pth = path to your model
%   maxclust = maximal number of clusters to find
%   method = clustering method
%   imgtype = 'spmT' or 'con' or 'beta'

function brainspmtConsistency(pth, maxclust, method, imgtype)
if nargin < 1
    pth = pwd;
end
if nargin < 2 || isempty(maxclust)
    maxclust = 2;
end
if nargin < 3 || isempty(method)
    method = 'complete';
end
if nargin < 4 || isempty(imgtype)
    imgtype = 'spmT';
end

colours = distinguishable_colors(maxclust);

D = dir(pth);

for d = length(D):-1:1
    if ~D(d).isdir || strcmp(D(d).name, 'diagnostics') || ...
            strcmp(D(d).name, '..') || strcmp(D(d).name, '.')
        D(d) = [];
    else
        % In case we have a folder within a folder...
        DD = dir(fullfile(pth, D(d).name));
        if length(DD) == 3 && DD(3).isdir
            D(d).name = fullfile(D(d).name, DD(3).name);
        end
    end
end

tD = dir(fullfile(pth, D(1).name, sprintf('%s_*.img', imgtype)));
if isempty(tD)
    tD = dir(fullfile(pth, D(1).name, sprintf('%s_*.nii',imgtype)));
end

for c = 1:length(tD)
    for d = 1:length(D)
        V = spm_vol(fullfile(pth, D(d).name, tD(c).name));
        Y = spm_read_vols(V);
        Y = Y(:);
        if d == 1
            data = nan(length(Y), length(D));
        end
        data(:,d) = Y;
    end
    M = all(isfinite(data) & data ~= 0,2);
    
    data = data(M, :);
    C = squareform(1-corrcoef(data));
    E = pdist(data');
    
    % Cluster participants' data...
    eZ = linkage(E, method);
    eT = cluster(eZ, 'maxclust', maxclust);
    eS = mdscale(E,2);
    
    cZ = linkage(C, method);
    cT = cluster(cZ, 'maxclust', maxclust);
    cS = mdscale(C,2);
    
    %% Plot everything...
    h.Fig = figure;
    subplot(2,3,1)
    imagesc(squareform(E));
    axis equal tight
    title(sprintf('Euclid (mn:%0.2f, sd:%0.2f)', mean(E), std(E)))
    
    subplot(2,3,4)
    imagesc(squareform(C));
    axis equal tight
    title(sprintf('1 - Correl. (mn:%0.2f, sd:%0.2f)', mean(C), std(C)))
    caxis([0 2])
    
    subplot(2,3,2);
    dendrogram(eZ);
    title('Dendrogram')
    
    subplot(2,3,5);
    dendrogram(cZ);
    title('Dendrogram')
    
    subplot(2,3,3);
    hold on
    for m = 1:maxclust
        scatter(eS(eT==m,1), eS(eT==m,2), [], colours(m,:));
    end
    title('Clusters')
    
    subplot(2,3,6);
    hold on
    for m = 1:maxclust
        scatter(cS(cT==m,1), cS(cT==m,2), [], colours(m,:));
    end
    title('Clusters')
    
    %% Find which voxels contribute most to similarity across participants
    R = 0;
    for x = 1:length(D)
        for y = 1:length(D)
            if y ~= x
                % This will only work with T-values... (option?)
                r = abs(data(:,y) - data(:,x));
                R = R + r;
                
                %{
                p = polyfit(data(:,x), data(:,y), 1);
                
                % Fitted y...
                f = p(1) .* data(:,x) + p(2);
                % Absolute residual...
                r = abs(data(:,y) - f);
                R = R + r;
                %}
                
                % DEBUG
                %{
                figure
                plot(data(:,x), data(:,y), 'x')
                hold on
                plot(data(:,x), f, '-');
                %}
            end
        end
    end
    I = nan(V.dim);
    I(M) = R;
    
    % Write consitency image...
    V.fname = fullfile(pth, sprintf('consistency_%04d.nii', c));
    spm_write_vol(V, I);
    
    % DEBUG
    % figure
    % plotmatrix_memspa(data);
end