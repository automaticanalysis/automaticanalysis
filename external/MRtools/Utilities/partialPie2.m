function [parts h] = partialPie2(y,x,labels)

if ~iscell(x);
    error('x should be cell array of sets');
end
y = zscore(y);

X = [];
for ii = 1:numel(x)
    x{ii} = zscore(x{ii});
    X = [X x{ii}];
end


TSS  = SumOfSquares(y);
Rall = SumOfSquares(X*(pinv(X)*y))/TSS;
res1 = y-(X*(pinv(X)*y));
b = (pinv(X)*y);
mm = SumOfSquares(X.*repmat(b',size(X,1),1));
SS = []; pSS = []; other = [];
for ii = 1:size(x,2);
    in = ii;
    out = setdiff(1:size(x,2),ii);
    
    X = [];
    for jj = out
        X = [X x{jj}];
    end
    
    R(ii) = SumOfSquares(X*(pinv(X)*y))/TSS;
    pR(ii) = Rall-R(ii);
end

%%
if sum(pR)>Rall
    warning('Suppression is in effect!');
    parts = [1-sum(pR) pR];
    extra = sum(pR)-Rall;
else
    parts = [1-Rall pR Rall-sum(pR)];
end


clf;
if nargin < 3
    h = pie(parts);
else
    leg = [];
    subplot(1,10,[1 2 3 4 6 7 ]);
    h = pie(parts); hold on;
    title(labels{1},'fontsize',18);
    leg{1} = ['Unexplained Variance = ' sprintf('%0.1f',parts(1)*100) '%'];
    
    
    for ii = 2:(numel(labels));
        leg{ii} = ['Unique to ' labels{ii} ' = ' sprintf('%0.1f',parts(ii)*100) '%'];
    end
   
    if ~exist('extra','var');
        hh = h(1:2:(numel(leg)*2));
        %leg{end+1} = ['Shared Variance with eYO = ' sprintf('%0.1f',parts(end)*100) '%'];
        leg{end+1} = ['Shared Variance = ' sprintf('%0.1f',parts(end)*100) '%'];
        hh(end+1) = h(end-1);
    else
        hh = h(1:2:(numel(leg)*2));
        th = plot(0,0);
        leg{end+1} = ['Enhancement = ' sprintf('%0.1f',extra*100) '%'];
        hh = [hh th];
    end
    delete(h(2:2:end));
    legend(hh,leg,'fontsize',14,'location','East');
end

