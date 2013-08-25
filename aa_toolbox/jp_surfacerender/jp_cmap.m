function cmap = jp_cmap(c, nsteps)
%JP_CMAP Create custom single-color colormap.
%
%  JP_CMAP(C) creates a custom colormap using C, which is an 1x3 vector
%  indicating a Matlab color.
%
%  JP_CMAP(C,NSTEPS) creates the colormap with NSTEPS number of steps.
%  If not specified, 64 steps are used.
%
%  Example:
%
%    red = jp_cmap([1 0 0]);

% Jonathan Peelle
% University of Pennsylvania

if nargin < 2 || isempty(nsteps)
  nsteps = 64;
end

if size(c,1) > 1 || size(c,2) ~= 3
  error('C must be a 1 by 3 vector indicating a Matlab color.');
end

cmap = repmat(c, nsteps, 1);




