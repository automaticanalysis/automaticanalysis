function ScatterPlot(X,Y,Labs,AxLabs)
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

h = [];
Leg = [];
pp = [];
for ii = 1:numel(X);
    if numel(X)>2
        subplot(1,numel(X),ii);
    end
    x = X{ii}; y = Y{ii};
    pp(ii) = plot(x,y,'o','markersize',8,'color','k','linewidth',2, 'markerfacecolor',cols(ii,:)); hold on;
    
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
    
    % xx = min(x):(max(x)-min(x))/100:max(x);
    color = cols(ii,:); color = (color+[.9 .9 .9])/2;

    ny = b(1)+(b(2)*xx);
    h(end+1) = plot(xx,ny,   '-','linewidth',4.0,'Color',color);
    h(end+1) = plot(xx,ny+ci,'-','linewidth',2.0,'Color',color);
    h(end+1) = plot(xx,ny-ci,'-','linewidth',2.0,'Color',color);

    r = sqrt(sum((pred-mean(pred)).^2)./sum((y-mean(y)).^2))*sign(b(2));
    Leg{end+1} = [Labs{ii} ': r = ' sprintf('%1.4f',r)];
    axis tight; axis square;
    xlabel(AxLabs{1});
    ylabel(AxLabs{2});
    
    if numel(X)>2
        title(Labs{ii})
    end
end
legend(pp,Leg,'location','best','fontsize',18);

for ii = 1:numel(h);
    uistack(h(ii),'bottom');
end
