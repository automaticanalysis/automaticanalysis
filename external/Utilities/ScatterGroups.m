function [h X] = ScatterGroups(dat,Groups,Labs)

tt = unique(Groups);
for ii = tt(:)'
    ind = find(Groups==ii);
    %plot(ii,dat(ind),'bd','markersize',10); hold on;
    
    mm(ii) = mean(dat(ind));
    sd(ii) = std(dat(ind));
    se(ii) = (std(dat(ind))./sqrt(numel(ind)))*1.96;
end

% errorbar(1:length(mm),mm,2*sd,'.k','linewidth',1.5)
eh = errorbar((1:length(mm)),mm,se,'o','color',[.65 .65 .65],'linewidth',4,'markersize',10,'markerfacecolor',[.65 .65 .65]); hold on;
errorbar_tick(eh,10); hold on

h = [];
for ii = tt(:)'
    ind = find(Groups==ii);
    X{ii} = ii+(.2*(rand(numel(ind),1)-.5));
    h(ii) = plot(X{ii},dat(ind),'b.','markersize',20); hold on;
end

bp = [];
if exist('boxplot.m')>0
%     [junk user] = UserTime;  
%     if strcmp(user,'aschultz');
%         bp = boxplot(dat,Groups,'notch','on');
%     else
        bp = boxplot(dat,Groups);
%     end
    set(findobj(gcf,'color', 'k'),'linewidth',3)
    set(findobj(gcf,'color', 'r'),'linewidth',3,'color','k')
    set(findobj(gcf,'color', 'b'),'linewidth',3,'color','k')
end

if numel(h)==1
    cols = [0 0 1];
else
    cols = zeros(numel(h),3);
    cols(:,3) = 1:-1/(numel(h)-1):0;
    cols(:,1) = 0:1/(numel(h)-1):1;
end
    
for ii = 1:numel(h)
    set(h(ii),'Color',cols(ii,:));
end

set(gca,'XTick', 1:length(Labs), 'XTickLabel',Labs,'FontSize',14,'FontWeight','bold');

uistack(eh,'bottom');
if ~isempty(bp)
    uistack(bp,'bottom');
end


ax = axis;
ax(1:2) = [.5 numel(tt)+.5];
axis(ax);

plot(ax(1:2),[0 0],'m--','linewidth',2);

function errorbar_tick(h,w,xtype)
%ERRORBAR_TICK Adjust the width of errorbars
%   ERRORBAR_TICK(H) adjust the width of error bars with handle H.
%      Error bars width is given as a ratio of X axis length (1/80).
%   ERRORBAR_TICK(H,W) adjust the width of error bars with handle H.
%      The input W is given as a ratio of X axis length (1/W). The result 
%      is independent of the x-axis units. A ratio between 20 and 80 is usually fine.
%   ERRORBAR_TICK(H,W,'UNITS') adjust the width of error bars with handle H.
%      The input W is given in the units of the current x-axis.
%
%   See also ERRORBAR
%

% Author: Arnaud Laurent
% Creation : Jan 29th 2009
% MATLAB version: R2007a
%
% Notes: This function was created from a post on the french forum :
% http://www.developpez.net/forums/f148/environnements-developpement/matlab/
% Author : Jerome Briot (Dut) 
%   http://www.mathworks.com/matlabcentral/newsreader/author/94805
%   http://www.developpez.net/forums/u125006/dut/
% It was further modified by Arnaud Laurent and Jerome Briot.

% Check numbers of arguments
error(nargchk(1,3,nargin))

% Check for the use of V6 flag ( even if it is depreciated ;) )
flagtype = get(h,'type');

% Check number of arguments and provide missing values
if nargin==1
	w = 80;
end

if nargin<3
   xtype = 'ratio';
end

% Calculate width of error bars
if ~strcmpi(xtype,'units')
    dx = diff(get(gca,'XLim'));	% Retrieve x limits from current axis
    w = dx/w;                   % Errorbar width
end

% Plot error bars
if strcmpi(flagtype,'hggroup') % ERRORBAR(...)
    
    hh=get(h,'children');		% Retrieve info from errorbar plot
    x = get(hh(2),'xdata');		% Get xdata from errorbar plot
    
    x(4:9:end) = x(1:9:end)-w/2;	% Change xdata with respect to ratio
    x(7:9:end) = x(1:9:end)-w/2;
    x(5:9:end) = x(1:9:end)+w/2;
    x(8:9:end) = x(1:9:end)+w/2;

    set(hh(2),'xdata',x(:))	% Change error bars on the figure

else  % ERRORBAR('V6',...)
    
    x = get(h(1),'xdata');		% Get xdata from errorbar plot
    
    x(4:9:end) = x(1:9:end)-w/2;	% Change xdata with respect to the chosen ratio
    x(7:9:end) = x(1:9:end)-w/2;
    x(5:9:end) = x(1:9:end)+w/2;
    x(8:9:end) = x(1:9:end)+w/2;

    set(h(1),'xdata',x(:))	% Change error bars on the figure
    
end
end
end
% return
% %%
% dat = rand(200,1);
% groups = ones(200,1);
% groups(51:100)=2; groups(101:150) = 3; groups(151:200)=4;
% 
% figure(20); clf;
% ll = {'Group1' 'Group2' 'Group3' 'Group4'};
% clf; ScatterGroups(dat,groups,ll)