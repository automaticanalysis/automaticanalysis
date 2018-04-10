function hands = ScatterPlot(X,Y,Labs,AxLabs,ind)
%%% Written by Aaron Schultz (aschultz@martinos.org)
%%%
%%% Copyright (C) 2012,  Aaron P. Schultz
%%%
%%% Supported in part by the NIH funded Harvard Aging Brain Study (P01AG036694) and NIH R01-AG027435 
%%%
%%% This program is free software: you can redistribute it and/or modify
%%% it under the terms of the GNU General Public License as published by
%%% the Free Software Foundation, either version 3 of the License, or
%%% any later version.
%%% 
%%% This program is distributed in the hope that it will be useful,
%%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%% GNU General Public License for more details.

%%% ID the point CallBack
% command = [ ...
% 'if strcmpi(get(gco,''type''),''line'');' ...
% '    p = get(gca,''CurrentPoint'');' ...
% '    ps = [get(gco,''XData'')'' get(gco,''YData'')''];' ...
% '    dist = sqrt(sum((ps - repmat(p(1,1:2),size(ps,1),1)).^2,2));' ...
% '    disp(find(dist == min(dist)));' ...
% '    ind = find(dist ~= min(dist));' ...
% '    UD = get(gcf,''UserData'');' ...
% '    if(size(UD.X,2)<1); return; end;'...
% '    cla; reset(gcf); ScatterPlot(UD.X,UD.Y,UD.Labs,UD.AxLabs,ind);' ...
% 'end;' ...
% ];
% % set(h,'ButtonDownFcn',command);

% if nargin<5 && numel(X)==1
%     ind = 1:size(X{1},1);
% end

hands = [];

if ~iscell(X)
    X = {X};
end
if ~iscell(Y)
    Y = {Y};
end

if nargin<3 || isempty(Labs)
    for ii = 1:numel(X);
        Labs{ii} = ['C' num2str(ii)];
    end
end
   
if nargin<4 || isempty(AxLabs)
    AxLabs = {'V1' 'V2'};
end

if numel(X)==1
    cols = [0 0 1];
else
    cols = zeros(numel(X),3);
    cols(:,3) = 1:-1/(numel(X)-1):0;
    cols(:,1) = 0:1/(numel(X)-1):1;
end

xrange = [];
for ii = 1:numel(X);
    xrange(ii,1:2) = [min(X{ii}) max(X{ii})];
end
xrange = [min(xrange(:,1)) max(xrange(:,2))];
% keyboard;
h = [];
Leg = [];
pp = [];
for ii = 1:numel(X);
    if numel(X)>2
        subplot(1,numel(X),ii);
    end
    x = X{ii}; y = Y{ii};
    pp(ii) = plot(x,y,'o','markersize',8,'color','k','linewidth',2, 'markerfacecolor',cols(ii,:)); hold on;
    hands.pp(ii) = pp(ii);
    %set(pp(ii),'ButtonDownFcn',command);
    set(gca,'FontSize',18,'fontweight','bold')
    
    b = pinv([ones(numel(x),1) x])*y;
    pred = [ones(numel(x),1) x]*b;
    res = y-pred;
    
    %xx = unique(x);
    xx = xrange(1):diff(xrange)/100:xrange(2);
    
    tcrit = icdf('t',1-.025,numel(x)-2);
    rse = sqrt(sum(res.^2)/(numel(x)-2));
    
    ci = sqrt(1*((1/numel(x))+(((xx-mean(x)).^2)./sum((x-mean(x)).^2))));
    ci = tcrit*rse*ci;
    
    disp([tcrit rse tcrit*rse]);
    
    % xx = min(x):(max(x)-min(x))/100:max(x);
    color = cols(ii,:); color = (color+[.9 .9 .9])/2;

    ny = b(1)+(b(2)*xx);
    h(end+1) = plot(xx,ny,   '-','linewidth',4.0,'Color',color);
    h(end+1) = plot(xx,ny+ci,'-','linewidth',2.0,'Color',color);
    h(end+1) = plot(xx,ny-ci,'-','linewidth',2.0,'Color',color);
    hands.fit = h(end-2:end);
    
    r = sqrt(sum((pred-mean(pred)).^2)./sum((y-mean(y)).^2))*sign(b(2));
    Leg{end+1} = [Labs{ii} ': r = ' sprintf('%1.4f',r)];
    axis square;
    xlabel(AxLabs{1});
    ylabel(AxLabs{2});
    
    if numel(X)>2
        title(Labs{ii})
    end
end
hands.leg = legend(h(1:3:end),Leg,'location','best','fontsize',18);

for ii = 1:numel(h);
    uistack(h(ii),'bottom');
end

UD.Labs = Labs;
UD.AxLabs = AxLabs;
UD.X = x;
UD.Y = y;
set(gcf,'UserData',UD);
