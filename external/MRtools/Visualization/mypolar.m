function h = mypolar(y,labels,sets,lims,highlight)
%  % y is the matrix of values you want to plot, I?m just making up data here for 9 measurements and 3 groups.
% y = rand(9,3)*20;  
% 
% % Set the labels for the measurements, should be in the same order as the rows in y
% labels = {'DMN' 'DAN' 'VAN' 'FPCN' 'Lang' 'Prec' 'pVis' 'xVis' 'MOT'};
% 
% % Set the labels for the groups, corresponding to the 3 columns in y
% sets = {'Group1' 'Group2' 'Group3'};
% % Create the plot
% figure(10); clf;
% h = mypolar(y,labels,sets,[0 30]);

% keyboard;
%%% may want to set this up with configureable range parameter
oy = y;
% keyboard;
if nargin<4 || isempty(lims);
    lims = [min(y(:)) max(y(:))];
end
y = (y-lims(1))/range(lims);

% x = (0 : (2*pi)/(size(y,1)) : 2*pi)';
x = linspace(0,2*pi,size(y,1)+1)';

figure(gcf); clf;
% wh = reshape(1:100,10,10);
% wh = wh(:,1:8);
% subplot(10,10,wh(:));

xx = (0 : (2*pi)/100 : 2*pi)';
lev = [.2 .4 .6 .8 1];
for ii = 1:numel(lev)
    h.radials(ii) = plot(sin(xx)*lev(ii),cos(xx)*lev(ii),'color',[.75 .75 .75],'linewidth',.5);
    hold on;
end

if nargin==5 && ~isempty(highlight)
        sc = (highlight-lims(1))/range(lims);
        h.radials(end+1) = plot(sin(xx)*sc,cos(xx)*sc,'color',[0 0 0],'linewidth',2);
end

axis([-1 1 -1 1]);
axis square;
axis off;
set(gcf,'color','w');

%%%%
for ii = 1:numel(x)
    spokes = [0 0; 0 1]*[cos(x(ii)) -sin(x(ii)); sin(x(ii)) cos(x(ii))];
    spokes2 = [0 0; 0 1.1]*[cos(x(ii)) -sin(x(ii)); sin(x(ii)) cos(x(ii))];
    
    if ii<numel(x)
        h.spokes(ii) = plot(spokes(:,1),spokes(:,2),'k-','color',[.75 .75 .75]);
        h.text(ii) = text(spokes2(2,1), spokes2(2,2), labels{ii},'Color','k','fontsize',18,'fontname','Helvetica');
        
        
        t1 = get(h.text(ii),'Extent');
        t2 = get(h.text(ii),'position');
        t2(1) = t2(1) - (t1(3)/2);
        set(h.text(ii),'position',t2);
        
%         set(th,{'Rotation'},num2cell(  ( (0)+tan(spokes(2,1)/spokes(2,2)))/(2*pi) *360  ));
%         set(th,{'Rotation'},num2cell(  ( (0) + tan(x(ii)))/(2*pi) *360  ))
        
%     
%     s.t(~ipos) = text(x(~ipos)*1.1, y(~ipos)*1.1, lbls(~ipos),'Color','w','fontsize',18,'fontname','Helvetica');
%     set(s.t(~ipos),{'Rotation'}, num2cell(theta(~ipos)'/tau*360 - 180),'Horiz','right')
    end
end
%%%%
tx = [0 -1; 0 1]*[cos(((2*pi)/360)*90) -sin(((2*pi)/360)*90); sin(((2*pi)/360)*90) cos(((2*pi)/360)*90)];
h.scale = plot(tx(:,1),tx(:,2));

np = 7; % should always be an odd number
x1 = linspace(tx(1,1),tx(2,1),np);
y1 = linspace(tx(1,2),tx(2,2),np);

v = linspace(lims(1),lims(2),ceil(np/2));
h.scale = plot(x1,y1,'+-','color',[0 0 0],'linewidth',3,'markersize',10);
ord = [[ceil(np/2):-1:1]' [ceil(np/2):1:np]'];
h.scaleText = [];
for ii = 1:size(ord,1)
    
    h.scaleText(end+1) = text(x1(ord(ii,1)), y1(ord(ii,1))+.05, sprintf('%0.3f',v(ii)),'Color','k','fontsize',18,'fontname','Helvetica');
    t1 = get(h.scaleText(end),'Extent');
    t2 = get(h.scaleText(end),'position');
    t2(1) = t2(1) - (t1(3)/2);
    set(h.scaleText(end),'position',t2);
    
    h.scaleText(end+1) = text(x1(ord(ii,2)), y1(ord(ii,2))+.05, sprintf('%0.3f',v(ii)),'Color','k','fontsize',18,'fontname','Helvetica');
    t1 = get(h.scaleText(end),'Extent');
    t2 = get(h.scaleText(end),'position');
    t2(1) = t2(1) - (t1(3)/2);
    set(h.scaleText(end),'position',t2);
end
%%%%
cols = colmap('cbq_set2',9);
% keyboard;
for ii = 1:size(y,2)
    a = sin(x(1:end-1)).*y(:,ii);
    b = cos(x(1:end-1)).*y(:,ii);

    a(end+1) = a(1);
    b(end+1) = b(1);
    
    h.lines(ii) = plot(a,b,'.-','linewidth',3,'markersize',40,'color',cols(mod(ii-1,9)+1,:),'Visible','off');
    h.fills(ii) = fill(a,b,cols(mod(ii-1,9)+1,:),'FaceAlpha',.33,'EdgeColor',cols(mod(ii-1,9)+1,:),'Marker', '.', 'MarkerSize',40,'linewidth',3);
    hold on;
end

% v = abs(max(y(:))*1.1);

h.leg = legend(h.lines,sets,'location', 'Best','fontsize',24);

uistack(h.scale,'top'); uistack(h.scaleText,'top');
shg
