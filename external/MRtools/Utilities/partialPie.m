function [parts h] = partialPie(y,x,labels)

y = zscore(y);
x = zscore(x);
C = ones(size(x,1),1);

%%
TSS  = SumOfSquares(y);
Rall = SumOfSquares([C x]*(pinv([C x])*y))/TSS;
res1 = y-([C x]*(pinv([C x])*y));
b = (pinv([C x])*y);
mm = SumOfSquares(x.*repmat(b(2:end)',size(x,1),1));
SS = []; pSS = []; other = [];
for ii = 1:size(x,2);
    in = ii;
    out = setdiff(1:size(x,2),ii);
    X = [C x(:,out)];
    
    
    R(ii) = SumOfSquares(X*(pinv(X)*y))/TSS;
    pR(ii) = Rall-R(ii);
    R(ii) = corr(y,x(:,in)).^2;
    pR2(ii) = partialcorr(y,x(:,in),x(:,out)).^2;
    
end
[pR; pR2; R]
[Rall sum(pR) sum(pR2) sum(R)]


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
%     subplot(1,10,[1 2 3 4 6 7 ]);
    h = pie(parts); hold on;
    title(labels{1},'fontsize',18);
    leg{1} = ['Unexplained Variance = ' sprintf('%0.1f',parts(1)*100) '%'];
    
    
    for ii = 2:(numel(labels));
        leg{ii} = ['Unique to ' labels{ii} ' = ' sprintf('%0.1f',parts(ii)*100) '%'];
    end
    %keyboard 
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
    legend(hh,leg,'fontsize',14,'location','EastOutside');
end

%% Need to figure out how to make this work.

% a = randn(100,1);
% b = randn(100,1);
% c = randn(100,1);
% d = randn(100,1);
% 
% x = a +   (randn(100,1)*.4);
% y = b +   (randn(100,1)*.4);
% z = a+b + (randn(100,1)*.4);
% 
% % figure(20); clf;
% % [parts h] = partialPie(x,[y z],{'x' 'y' 'z'});
% 
% clc
% TSS = SumOfSquares(y);
% v1 = SumOfSquares(PredictedData(x,y))/TSS;
% v2 = SumOfSquares(PredictedData(x,z))/TSS;
% v3 = SumOfSquares(PredictedData(x,[y z]))/TSS;
% 
% [v3 v1 v2 v3-v1 v3-v2]
% 
% 
% % ([partialcorr(x,y,z) partialcorr(x,z,y) corr(x,z) corr(x,y)])




