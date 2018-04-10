function [h, pp, l] = ScatterPlot2(X,Y,Groups,AxLabs,sep)
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
% 'end;' ...
% ];
% set(gcf,'ButtonDownFcn',command);

if nargin<4 || isempty(AxLabs)
    AxLabs = {'V1' 'V2'};
end

if nargin<5
    sep=false;
end


if isobject(Groups)
    Groups = cellstr(Groups);
end

if isnumeric(Groups);
    Groups = cellstr(num2str(Groups(:)));
end

list = unique(Groups,'stable');

Labs = list;

h = [];
Leg = [];
pp = [];
  
  
if numel(list)==1
    cols = [0 0 1];
else
%     cols = zeros(numel(list),3);
%     cols(:,3) = 1:-1/(numel(list)-1):0;
%     cols(:,1) = 0:1/(numel(list)-1):1;
    
    cols = colormap;
    wh =linspace(1,size(cols,1),numel(list));
    cols = cols(wh,:);
end
% keyboard;
xrange = [min(X) max(X)];
xx = xrange(1):diff(xrange)/100:xrange(2);

for ii = 1:numel(list);
    ind = strmatch(list{ii},Groups,'exact');
     
    if numel(list)>2 || sep
        subplot(1,numel(list),ii);
%         subplot(floor(numel(list)),ceil(sqrt(numel(list))),ii);
%         subplot(2,2,ii);
    end
    x = X(ind); y = Y(ind);
    if sep 
        xrange = [min(x) max(x)];
        xx = xrange(1):diff(xrange)/100:xrange(2);
    end
    pp(ii) = plot(x,y,'o','markersize',8,'color','k','linewidth',2, 'markerfacecolor',cols(ii,:)); hold on;
    %set(pp(ii),'ButtonDownFcn',command);
    set(gca,'FontSize',18,'fontweight','bold')
    
%     xx = xrange(1):diff(xrange)/100:xrange(2);
    xx = linspace(min(x),max(x),20);
    
    b = pinv([ones(numel(x),1) x])*y;
    pred = [ones(numel(x),1) x]*b;
    
    res = y-pred;   
    
    tcrit = icdf('t',1-.025,numel(x)-2);
    rse = sqrt(sum(res.^2)/(numel(x)-2));
%     disp([rse tcrit]);
    
    ci = sqrt(1*((1/numel(x))+(((xx-mean(x)).^2)./sum((x-mean(x)).^2))));
    ci = tcrit*rse*ci;
    
    color = cols(ii,:); color = (color+[.9 .9 .9])/2;
    ny = b(1)+(b(2)*xx);
    h(end+1) = plot(xx,ny,   '-','linewidth',4.0,'Color',color);
    h(end+1) = plot(xx,ny+ci,'-','linewidth',2.0,'Color',color);
    h(end+1) = plot(xx,ny-ci,'-','linewidth',2.0,'Color',color);

%     Z = smooth(x,y,.75,'rloess');
%     [~,ord] = sort(x);
%     color = cols(ii,:); color = (color+[.9 .9 .9])/2;
%     h(end+1) = plot(x(ord),Z(ord),   '-','linewidth',4.0,'Color',color);

    r = sqrt(sum((pred-mean(pred)).^2)./sum((y-mean(y)).^2))*sign(b(2));
    Leg{end+1} = [Labs{ii} ': r = ' sprintf('%1.3f',r)];
    axis square;
    xlabel(AxLabs{1});
    ylabel(AxLabs{2});
    
    if numel(list)>2 || sep
        title(Labs{ii})
    end
end
l = legend(h(1:3:end),Leg,'location','best','fontsize',18);

for ii = 1:numel(h);
    uistack(h(ii),'bottom');
end
