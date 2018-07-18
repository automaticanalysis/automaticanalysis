function h = SchemaBall(r, lbls, clims,colmap,opt,sets)
% SCHEMABALL Plots correlation matrix as a schemaball
%
%   SCHEMABALL(R) R is a square numeric symetric matrix].
%
%
%   SCHEMABALL(..., LBLS, CLIME, COLORMAP ,OPT, NCOLOR) Plot schemaball with optional 
%                                         arguments (accepts empty args).
%                                               
%       - LBLS      Plot a schemaball with custom labels at each node.
%                   LBLS is either a cellstring of length M, where
%                   M = size(r,1), or a M by N char array, where each
%                   row is a label.
%
%       - CLIMS     A 2 value vector denoting the minimum and maximum for
%                   the colorscale.
%      
%       - COLMAP    A string specify the desired colormap (see colmap.m for
%                   options)
%
%       - OPT       1 = Use quadratic bezier curves for all curves. 2 =
%                   Modify the exponent of the bezier curves to reflect
%                   proximity on the circle such that points cluser on the
%                   circle are connected by lower amplitude curves.
%
%       - SETS      An option parameter that will group items in labls into
%                   contiguous sets, and create an outer ring to help
%                   distinguish nodes, and color the labels according to
%                   the outer ring (will be done using the jet colormap)
%
%
%   H = SCHEMABALL(...) Returns a structure with handles to the graphic objects
%
%       h.l     handles to the curves (line objects), one per color shade. 
%               If no curves fall into a color shade that handle will be NaN.
%       h.s     handle  to the nodes (scattergroup object)
%       h.t     handles to the node text labels (text objects)
%       h.o     handles for the outer ring (only applicable when SETS is specified.
%       h.cb    handle for the colorbar.

% Examples
%
%   % Base demo
%   
%
%   % Supply your own correlation matrix (only lower off-diagonal triangular part is considered)
%   x = 2*(rand(10,10)-.5); 
%   SchemaBall(x,[],[-1 1],'jet',2,{1:3, 4:6 7:10})
%   
%   Now threshold x and change the color limits:
%
%   x(x<.5) = NaN;
%   SchemaBall(x,[],[.5 1],'jet',2,{1:3, 4:6 7:10})
%
%   This script is a modified version of schemaball.m written by Oleg
%   Komarov (oleg.komarov@hotmail.it) from 15 jun 2013, which was a modification of code
%   originally posted by Gunther Struyf.
%
%   Oleg Komarov's code is available through the MATLAB FEX
%   Gunther Struyf's code is available at https://github.com/GuntherStruyf/matlab-tools/blob/master/schemaball.m
%   A discussion thresh of these functions is available at http://www.stackoverflow.com/questions/17038377/how-to-visualize-correlation-matrix-as-a-schemaball-in-matlab/17111675
% 
%   
% 
%   This version was created by by Aaron Schultz (aschultz@martinos.org)
%   April 2-14
%
%   Copyright (C) 2011,  Aaron P. Schultz
%
%   Supported in part by the NIH funded Harvard Aging Brain Study (P01AG036694) and NIH R01-AG027435 
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.

t      = (0: 0.025 :1)';

% Some defaults
sz = size(r);
if nargin < 2 || isempty(lbls);     lbls   = cellstr(reshape(sprintf('%-4d',1:sz(1)),4,sz(1))');  end

% Lbls
if (~ischar(lbls) || size(lbls,1) ~= sz(1)) && (~iscellstr(lbls) || ~isvector(lbls) || length(lbls) ~= sz(1))
    error('schemaball:validLbls','LBLS should either be an M by N char array or a cellstring of length M, where M is size(R,1).')
end
if ischar(lbls)
    lbls = cellstr(lbls);
end

%%% Set up the colomap and color limits
% r = tril(r,-1);
% r(r==0)=NaN;
i1 = find(~isnan(r));
cols = nan(numel(r),3);
cc = nan(numel(r),1);
[cols1 cm cc1] = cmap(r(i1), clims, colmap);
cols(i1,:) = cols1;
cc(i1) = cc1;

%%% Engine
% Create figure
% figure('renderer','zbuffer','visible','off') % ,'Renderer','OpenGL'
clf;
axes('NextPlot','add')
% set(gcf,'Visible','off');
% Index only low triangular matrix without main diag
tf        = tril(true(sz),-1);

%%%
% Retrieve pairings of nodes
[row,col] = find(tf);

% Use tau http://tauday.com/tau-manifesto
tau   = 2*pi;
% Positions of nodes on the circle starting from (0,-1), useful later for label orientation
step  = tau/sz(1);
theta = -.25*tau : step : .75*tau - step;
% Get cartesian x-y coordinates of the nodes
x     = cos(theta);
y     = sin(theta);
%%% PLOT BEZIER CURVES 
if opt == 1
    %%%  use a constant exponent
    t2  = [1-t, t].^2;
    N2 = unique(cc1);
    s.l = NaN(numel(N2),1);
    % LOOP per color bucket
    for c = 1:numel(N2)
        
        pos = i1(find(cc1==N2(c)));
        [row,col] = ind2sub(size(r),pos);
        
        Bx     = [t2*[x(col); x(row)]; NaN(1,numel(pos))]; % NaN(1,numel(pos))
        By     = [t2*[y(col); y(row)]; NaN(1,numel(pos))]; % NaN(1,numel(pos))
        s.l(c) = plot(Bx(:),By(:),'Color',cols(pos(1),:),'LineWidth',2); % round(sqrt(cc(pos))/4)   ,'linesmoothing','on'
    end
elseif opt == 2
    %%% Plot each curve separately using different exponent for each curve
    %s.l = NaN(sum(~isnan(cc)),1);
    tmp = tril(r,-1); tmp(tmp==0)=NaN; 
    ind = find(~isnan(tmp));
    [trash i] = sort(r(ind));
    ind = ind(i);
    s.l = nan(numel(ind),1);
    s.li = i;
    cc = 0;
    for c = ind(:)'
        cc = cc+1;
        [rw,cl] = ind2sub(sz,c);
        dist = (sqrt(((x(rw)-x(cl))^2 + (y(rw)-y(cl))^2))/2)*2;
        
        t2  = [1-t, t].^(1+dist);
        
        Bx     = [t2*[x(cl); x(rw)]; ]; 
        By     = [t2*[y(cl); y(rw)]; ]; 
        s.l(cc) = plot(Bx(:),By(:),'Color',cols(c,:),'LineWidth',2); % round(sqrt(cc(pos))/4) ,'linesmoothing','on'
%         set(s.l(end),'LineWidth',9/(256./cc(c)));
%         keyboard;
    end
end

%%%
% PLOT NODES
vals = sum(~isnan(r))';
% if any(r(:)<1)
%     vals = nanmean(r.^2)';
% else
%     vals = nanmean(r)';
% end

% i2 = find(~isnan(vals)); i2 = i2(:);
% kk = zeros(numel(vals),3);
% [k1 k2 k3] = cmap(vals(i2),[0 max(vals)],colmap);
% kk(i2,:) = k1;

[k1 k2 k3] = cmap(vals,[0 max(vals)],colmap);
kk = k1;

% Plot in brighter color those nodes which on average are more absolutely correlated
s.s = scatter(x,y,100,kk,'fill','MarkerEdgeColor','w','linewidth',2); %'MarkerEdgeColor',ecolor,'LineWidth',1
%%%
% PLACE TEXT LABELS such that you always read 'left to right'
ipos       = x > 0;
s.t        = zeros(sz(1),1);
s.t( ipos) = text(x( ipos)*1.1, y( ipos)*1.1, lbls( ipos),'Color','w','fontsize',18,'fontname','Helvetica');
set(s.t( ipos),{'Rotation'}, num2cell(theta(ipos)'/tau*360))
s.t(~ipos) = text(x(~ipos)*1.1, y(~ipos)*1.1, lbls(~ipos),'Color','w','fontsize',18,'fontname','Helvetica');
set(s.t(~ipos),{'Rotation'}, num2cell(theta(~ipos)'/tau*360 - 180),'Horiz','right')


set(gca,'position',[0 .1 1 .8]);
axis equal; axis off; shg

xtn        = cell2mat(get(s.t,'extent'));
post       = cell2mat(get(s.t,'pos'));
sg         = sign(post(:,2));
posfa      = cell2mat(get([gcf gca],'pos'));
% Calculate xlim and ylim in data units as x (y) position + extension along x (y)
ylims      = post(:,2) + xtn(:,4).*sg;
ylims      = [min(ylims), max(ylims)];
xlims      = post(:,1) + xtn(:,3).*sg;
xlims      = [min(xlims), max(xlims)];


val = max(abs([xlims ylims]));
% prevent resize
axis([-val val -val val]);

set(gca,'color','k');
set(gcf,'color','k');


sc = val*1.1;
%%%%
if nargin>5 && ~isempty(sets);
    
    tau   = 2*pi;
    step  = tau/sz(1); 
    theta = (-.25*tau) : step : .75*tau ;% -step;
    % Get cartesian x-y coordinates of the nodes
    x     = cos(theta);
    y     = sin(theta);
    
    cc = jet(numel(sets));
    hh = nan(numel(sets,1));

    for ii = 1:numel(sets)       
        ind = sets{ii};
        a = theta(ind(1)) -(step*.4);
        b = theta(ind(end)) + (step*.4);
        x = a:(b-a)/50:b;
        hh(ii) = plot(cos(x)*sc*1.1,sin(x)*sc*1.1,'-','color',cc(ii,:),'linewidth',15); %,'linesmoothing','on'
        set(s.t(sets{ii}),'color',cc(ii,:));
%         shg; pause
    end
    axis([-sc sc -sc sc]*1.2);
    s.o = hh;
end

%%%%
colormap(cm);
cb = colorbar;
h.cb = cb;
xx = clims(1):diff(clims)/4:clims(2);
labs = {};
for ii = 1:numel(xx);
    labs{ii} = sprintf('%0.3f',xx(ii));
end
set(cb,'YTick',0:1/(numel(xx)-1):1, 'YTickLabel',labs, 'fontsize',18,'fontname','Helvetica','color','w');
set(gca,'fontname','Helvetica');
%%%%
h = s;
shg

% set(gcf,'Visible','on');