% SCATTER2HEAT
% scatterMatrix = scatter2heat(X, Y, C, bins, transform)
% Plots a memory efficient scatterplot by creating a heatmap
% Adapted from mathematical.coffee at stackoverflow.com
%
% Inputs:
% X: x axis points of scatterplot
% Y: y axis points of scatterplot
% C: how to colour the points (if empty, the heatmap will contain number of
% points in that bin)
% bins: how many bins in each direction...
function scatterMatrix = scatter2heat(X, Y, C, bins, transform)

if nargin < 4
    % Group into bins of 0.05
    bins = 20;
end
if nargin < 5
    % Group into bins of 0.05
    transform = '';
end
if length(bins) == 1
    bins = [bins, bins];
end

% Create the bins
xbins = linspace(min(X), max(X), bins(1));
ybins = linspace(min(Y), max(Y), bins(2));

% Work out which bin each X and Y is in (idxx, idxy)
[nx, idxx] = histc(X, xbins);
[ny, idxy] = histc(Y, ybins);

% Calculate mean (or number of points) in each direction
if nargin < 3 || isempty(C);
    scatterMatrix = accumarray([idxy idxx], ones(size(X))', [], @sum);
else
    scatterMatrix = accumarray([idxy idxx], C', [], @mean);
end

if ~isempty(transform)
    eval(['scatterMatrix = ' transform '(scatterMatrix);']);
end

% PLOT
imagesc(xbins, ybins, scatterMatrix)
set(gca, 'YDir', 'normal') % flip Y axis back to normal
colorbar