function jp_spm8_surfacerender(img, cmap, cfg)
%JP_SPM8_SURFACERENDER render some data on cortical surfaces.
%
% JP_SPM8_SURFACERENDER(IMG) renders the data from IMG onto surfaces. If
% not specified or left empty, you will be prompted to select one.
%
% If IMG is a string, and not an image, it gets passed to SPM_SELECT as a
% filter.  For example JP_SPM8_SURFACERENDER('^spmT.*') would make it
% easier to select T images.
%
% JP_SPM8_SURFACERENDER(IMG, COLORMAP) lets you specify the colormap used
% (default 'hot').  You can also specify a single color as an RGB triplet
% (e.g., [.5 .4 0]).
%
% JP_SPM8_SURFACERENDER(IMG, COLORMAP, CFG) lets you specify additional
% options.
%
% For the colorscale, leaving cfg.colorscale empty just does the default,
% which scales from min to max in the data.  'Symmetrical' tries to make it
% so 0 is the middle of the color scale (so using a jet map, for example,
% green = 0).  Specifying [min max] lets you keep consistent scaling across
% several figures (instead of rescaling t maps arbitrarily); values below
% the min are not displayed.
%
% cfg = [];
% cfg.colorscale = [2 5]; % 'symmetrical' | [min max] | []
%
% cfg.plots = [1 2 3 4];  % specify a subset if you don't want all 4 views
%                         % 1=left, 2=right, 3=top, 4=bottom
%
%
% cfg.inflate = 5;        % # bigger number=more inflated. 0=don't inflate
% cfg.figname = 'my name'; % name the figure window
%
%
% There are also options relating to the figure size and placement of plots
% which allow you to customize the output:
%
% cfg.plot1pos = [X Y width height];
% cfg.plot2pos = [.42 .4 .35 .3] % if only using L and R later views
%
%
% Finally, there is a modified version of spm_mesh_project, which just
% displays the whole image, instead of > 0.  This is called
% spm_mesh_project_jp, and just needs to be in your path somewhere (e.g.,
% the same folder as JP_SPM8_SURFACERENDER).

% Jonathan Peelle
% University of Pennsylvania

% (just to remind myself what we need in the dat structure)
% dat      - a struct array of length 1 to 3
%            each element is a structure containing:
%            - XYZ - the x, y & z coordinates of the transformed SPM{.}
%                    values in units of voxels.
%            - t   - the SPM{.} values.
%            - mat - affine matrix mapping from XYZ voxels to MNI.
%            - dim - dimensions of volume from which XYZ is drawn.

if nargin < 1
  img = [];
end
  
  
if isempty(img) || ~exist(img)
  % If not an image, assume 'img' is a filter for the image
  if ischar(img) && ~isempty(img)
    img = spm_select(1, 'image', 'Select file with data to render', [], [], img);
  else    
    img = spm_select(1, 'image', 'Select file with data to render');
  end
end

if nargin < 2 || isempty(cmap)
  cmap = 'hot';
end

if nargin < 3
  cfg = [];
end

if ~isfield(cfg, 'plots') || isempty(cfg.plots)
  cfg.plots = [1:4];
end

if ~isfield(cfg, 'plot1pos') || isempty(cfg.plot1pos)
  cfg.plot1pos = [.05 .4 .35 .3];
end

if ~isfield(cfg, 'plot2pos') || isempty(cfg.plot2pos)
  cfg.plot2pos = [.6 .4 .35 .3];
end

if ~isfield(cfg, 'plot3pos') || isempty(cfg.plot3pos)
  cfg.plot3pos = [.37 .5 .25 .2];
end

if ~isfield(cfg, 'plot4pos') || isempty(cfg.plot4pos)
  cfg.plot4pos = [.37 .25 .25 .2];
end

if ~isfield(cfg, 'figname')
    cfg.figname = '';
end


if ~isfield(cfg, 'inflate') || isempty(cfg.inflate)
  cfg.inflate = 5;
end

if ~isfield(cfg, 'colorscale')
  cfg.colorscale = [];
end

if ~isfield(cfg, 'rendfile') || isempty(cfg.rendfile)
  cfg.rendfile = fullfile(spm('dir'), 'canonical', 'cortex_20484.surf.gii');
end


% get the surface file to render on
rend = export(gifti(cfg.rendfile),'patch');

% note what image we're rendering
fprintf('Rendering %s\n', img);

% set the colormap
if ischar(cmap)
    col = eval(sprintf('%s(256);', cmap));
else    
    col = jp_cmap(cmap, 256); % make a colormap that SPM is happy with
end
%col = eval(sprintf('%s(256);', cmap));

% get the data from the image
V = spm_vol(img);
[Y,XYZ] = spm_read_vols(V);
XYZvox = round(V.mat\[XYZ; ones(1,size(XYZ,2))]);

fprintf('Range of raw data between %g and %g.\n', min(min(min(Y))), max(max(max(Y))));

% if [min max] colorscale specified, trim data to only display > min value
if ~isempty(cfg.colorscale) && ~strcmp(cfg.colorscale, 'symmetrical')
  Ysiz = size(Y);
  YY = reshape(Y,1,numel(Y));
  YY(YY < cfg.colorscale(1)) = NaN;  
  Y = reshape(YY, Ysiz);
  clear YY Ysiz
end


dat(1) = struct( 'XYZ',  XYZvox,...
  't',    Y,...
  'mat',  V.mat,...
  'dim',  V.dim');



fprintf('Range of rendered data between %g and %g.\n', min(min(min(Y))), max(max(max(Y))));

% create a blank window to put these brains in
myfig = figure('color', 'w', 'position', [79 72 1034 850], 'name', cfg.figname);




%-------------------------------------------------------------
% Render the surfaces
%-------------------------------------------------------------

% surf_rend(dat, rend, colormap, figurename, which_brain, cfg)
% which_brain: 1,2,3,4 = L,R,Top,Bottom

if ismember(1, cfg.plots)
  hp1 = surf_rend(dat,rend,col,myfig, 1, cfg); % L
end

if ismember(2, cfg.plots)
  hp2 = surf_rend(dat,rend,col,myfig, 2, cfg); % R
end

if ismember(3, cfg.plots)
  hp3 = surf_rend(dat,rend,col,myfig, 3, cfg); % Top
end

if ismember(4, cfg.plots)
  hp4 = surf_rend(dat,rend,col,myfig, 4, cfg); % Bottom
end



%-------------------------------------------------------------
% Do some inflating
%-------------------------------------------------------------

if cfg.inflate > 0
  
  fprintf('Inflating...')
  
  % spm_mesh_inflate(handle_to_surface, how_many_steps, update_every_how_many)
  % (I like it about 5 steps inflated)
  
  if ismember(1, cfg.plots)
    spm_mesh_inflate(hp1,cfg.inflate,cfg.inflate);
  end
  
  if ismember(2, cfg.plots)
    spm_mesh_inflate(hp2,cfg.inflate,cfg.inflate);
  end
  
  if ismember(3, cfg.plots)
    spm_mesh_inflate(hp3,cfg.inflate,cfg.inflate);
  end
  
  if ismember(4, cfg.plots)
    spm_mesh_inflate(hp4,cfg.inflate,cfg.inflate);
  end
  
  fprintf('done.\n\n')
end

fprintf('\n      ''\n     ''\n  ____''___\n  |      |_\n  |      |_|\n  |______|    All done!\n\n'); % coffee




%-------------------------------------------------------------
% (below here are hacked versions of what was here to begin with)
%-------------------------------------------------------------


%==========================================================================
% function surf_rend(dat,rend,col)
%==========================================================================
function hp = surf_rend(dat,rend, col, myfig, plotnum, cfg)

%-Setup figure and axis
%--------------------------------------------------------------------------
Fgraph = myfig; %spm_figure('GetWin','Graphics');
%spm_results_ui('Clear',Fgraph);
rdr = get(Fgraph,'Renderer');
set(Fgraph,'Renderer','OpenGL');

% switch plotnum
%   case 1
%     ax = axes('position', [.05 .4 .35 .3], 'visible', 'off');
%     myview = [-90 0]; % left
%   case 2
%     % if only doing L and R hemisphere, space differently
%     if length(cfg.plots)==2 && sum(cfg.plots==[1 2])==2
%       fprintf('hiiii!!!!!!\n\n\n\n;');
%       ax = axes('position', [.42 .4 .35 .3], 'visible', 'off');
%     else
%       ax = axes('position', [.6 .4 .35 .3], 'visible', 'off');
%     end
%     myview = [90 0]; % right
%   case 3
%     ax = axes('position', [.37 .5 .25 .2], 'visible', 'off');
%     myview = [0 90]; % top
%   case 4
%     ax = axes('position', [.37 .25 .25 .2], 'visible', 'off');
%     myview = [-180 -90]; %bottom
%   otherwise
%     error('Only know about 4 types of plot positions.');
% end

switch plotnum
  case 1
    ax = axes('position', cfg.plot1pos, 'visible', 'off');
    myview = [-90 0]; % left
  case 2
    % if only doing L and R hemisphere, space differently
    if length(cfg.plots)==2 && sum(cfg.plots==[1 2])==2
      ax = axes('position', cfg.plot2pos, 'visible', 'off');
    else
      ax = axes('position', cfg.plot2pos, 'visible', 'off');
    end
    myview = [90 0]; % right
  case 3
    ax = axes('position', cfg.plot3pos, 'visible', 'off');
    myview = [0 90]; % top
  case 4
    ax = axes('position', cfg.plot4pos, 'visible', 'off');
    myview = [-180 -90]; %bottom
  otherwise
    error('Only know about 4 types of plot positions.');
end


%-Project data onto surface mesh
%--------------------------------------------------------------------------
v = spm_mesh_project_jp(rend, dat);  % jp modification - show both positive and negative

%-Compute mesh curvature texture
%--------------------------------------------------------------------------
curv = spm_mesh_curvature(rend) > 0;
curv = 0.5 * repmat(curv,1,3) + 0.3 * repmat(~curv,1,3);

%-Combine projected data and mesh curvature
%--------------------------------------------------------------------------
cdat = zeros(size(v,2),3);

% make equal around 0
maxv = max(abs(v)); % largest positive or negative value
minv = -1 * maxv;


if any(v(:))
  if strcmp(cfg.colorscale, 'symmetrical')
    cdat = squeeze(ind2rgb(floor( (v(:)-minv)/(2*maxv) * size(col,1)), col));
    fprintf('Color scale from -%g to %g.\n', maxv, maxv);
  elseif ~isempty(cfg.colorscale)
    cmin = cfg.colorscale(1);
    cmax = cfg.colorscale(2);
    
    cdat = squeeze(ind2rgb(floor( (v(:)-cmin)/(cmax) * size(col,1)), col));
    fprintf('Color scale from %g to %g.\n', cmin, cmax);    
  else
    if size(col,1)>3
      cdat = squeeze(ind2rgb(floor(v(:)/max(v(:))*size(col,1)),col));
    else
      m = max(v(:));
      for i=1:size(v,1)
        cdat = cdat + v(i,:)'/m * col(i,:);
      end
    end
  end
end  
 
% if any(v(:))
%     if size(col,1)>3
%         cdat = squeeze(ind2rgb(floor(v(:)/max(v(:))*size(col,1)),col));
%     else
%         m = max(v(:));
%         for i=1:size(v,1)
%             cdat = cdat + v(i,:)'/m * col(i,:);
%         end
%     end
% end

cdat = repmat(~any(v,1),3,1)' .* curv + repmat(any(v,1),3,1)' .* cdat;

%-Display the surface mesh with texture
%--------------------------------------------------------------------------
hp = patch(rend, 'Parent',ax,...
    'FaceVertexCData',cdat, ...
    'FaceColor', 'interp', ...
    'EdgeColor', 'none',...
    'FaceLighting', 'phong',...
    'SpecularStrength' ,0.7, 'AmbientStrength', 0.1,...
    'DiffuseStrength', 0.7, 'SpecularExponent', 10,...
    'DeleteFcn', {@mydeletefcn,Fgraph,rdr});

view(ax,myview);
axis(ax,'image');

l = camlight; set(l,'Parent',ax);
material(Fgraph,'dull');
setappdata(ax,'camlight',l);


%==========================================================================
function myinflate(obj,evd)
spm_mesh_inflate(getappdata(obj,'patch'),Inf,1);
axis(getappdata(obj,'axis'),'image');


%==========================================================================
function myview(obj,evd,varargin)
view(getappdata(get(obj,'parent'),'axis'),varargin{1});
axis(getappdata(get(obj,'parent'),'axis'),'image');
camlight(getappdata(getappdata(get(obj,'parent'),'axis'),'camlight'));

%==========================================================================
function mytransparency(obj,evd)
t = 1 - sscanf(get(obj,'Label'),'%d%%') / 100;
set(getappdata(get(obj,'parent'),'patch'),'FaceAlpha',t);
set(get(get(obj,'parent'),'children'),'Checked','off');
set(obj,'Checked','on');


%==========================================================================
function mypostcallback(obj,evd)
try, camlight(getappdata(evd.Axes,'camlight')); end

%==========================================================================
function mydeletefcn(obj,evd,varargin)
try, rotate3d(get(obj,'parent'),'off'); end
set(varargin{1},'Renderer',varargin{2});
