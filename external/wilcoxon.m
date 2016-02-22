function STATS=wilcoxon(varargin)
%This file execute the non parametric Wilcoxon test to evaluate the difference between paired (dependent) samples. 
%If the number of difference is less than 15, the algorithm calculate the exact ranks distribution; 
%else it uses a normal distribution approximation. 
%Now, the MatLab function SIGNRANK returns the same p-value. 
%Anyway, this Wilcoxon function gives a more detailed output (that is necessary for publications...)
%By itself Wilcoxon will run the demo
%
% Syntax: 	STATS=WILCOXON(X1,X2,PLTS)
%      
%     Inputs:
%           X1 and X2 - data vectors.
%           PLTS - Flag to set if you don't want (0) or want (1) view the plots
%     Outputs:
%           - W value and p-value when exact ranks distribution is used.
%           - W value, Z value, Standard deviation (Mean=0), p-value when normal distribution is used
%        If STATS nargout was specified the results will be stored in the STATS
%        struct.
%
%      Example: 
%
%         X1=[77 79 79 80 80 81 81 81 81 82 82 82 82 83 83 84 84 84 84 85 85 86 86 87 87];
% 
%         X2=[82 82 83 84 84 85 85 86 86 86 86 86 86 86 86 86 87 87 87 88 88 88 89 90 90];
%
%           Calling on Matlab the function: wilcoxon(X1,X2)
%
%           Answer is:
%
% WILCOXON TEST
% --------------------------------------------------------------------------------
% Sample size is good enough to use the normal distribution approximation
%  
% W         mW		sW          zW          p-value (2-tailed)
% --------------------------------------------------------------------------------
% 325  	     0		73.1608   4.4354	0.0000
% --------------------------------------------------------------------------------
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2006). Wilcoxon test: non parametric Wilcoxon test for paired samples.
% http://www.mathworks.com/matlabcentral/fileexchange/12702

%Input Error handling
args=cell(varargin);
nu=numel(args);
default.values = {[77 79 79 80 80 81 81 81 81 82 82 82 82 83 83 84 84 84 84 85 85 86 86 87 87],[82 82 83 84 84 85 85 86 86 86 86 86 86 86 86 86 87 87 87 88 88 88 89 90 90],0};
default.values(1:nu) = args;
[x1 x2 plts] = deal(default.values{:});
if nu==0
    plts=1;
elseif nu==1
    error('Warning: Almost two input data are required')
elseif nu>=2
    if ~isvector(x1) || ~isvector(x2)
        error('WILCOXON requires vector rather than matrix data.');
    end
    if (numel(x1) ~= numel(x2))
        error('Warning: WILCOXON requires the data vectors to have the same number of elements.');
    end
    if ~all(isfinite(x1)) || ~all(isnumeric(x1)) || ~all(isfinite(x2)) || ~all(isnumeric(x2))
        error('Warning: all X1 and X2 values must be numeric and finite')
    end
elseif nu==3
    if plts ~= 0 && plts ~= 1 %check if plts is 0 or 1
        error('Warning: PLTS must be 0 if you don''t want or 1 if you want to see plots.')
    end
else
    error('Warning: Max three input data are required')
end 
clear args default nu


tr=repmat('-',1,80);
disp('WILCOXON TEST')
disp(tr)
dff=sort(x2-x1); %difference between x1 and x2
dff(dff==0)=[]; %eliminate null variations
n=length(dff); %number of ranks
if length(x1)~=n %tell me if there are null variations
    fprintf('There are %d null variations that will be deleted\n',length(x1)-n)
end
if isempty(dff) %if all variations are null variations exit function
    disp('There are not variations. Wilcoxon test can''t be performed')
    return       
end

%Ranks of absolute value of samples differences with sign
[r,t]=tiedrank(abs(dff)); %ranks and ties
W=sum(r.*sign(dff)); %Wilcoxon statics (sum of ranks with sign)
pem=median(dff); %point estimation of median of differences
m=ceil(n/2); %location of the median
if mod(n,2)==0 %If the length of the series is even
    tmp=[dff(1:m) pem dff(m:end)]; %add the median in the middle
    dff=tmp; clear tmp 
    m=m+1;
end
%find how many differences far from the median we have to choose
C=cumsum(binopdf(0:1:n,n,0.5));
T=find(C<=0.025,1,'last')-1;
cintpem=dff([m-T m+T]); %construct the interval
clear C T m
fprintf('Median of differences (binomial estimator): %0.4f\n',pem)
fprintf('95%% Confidence interval: %0.4f  %0.4f\n',cintpem(1),cintpem(2))
disp('')
[I J]=ndgrid(dff,dff); d=triu(I+J)./2; %Walsh averages triangular matrix
ld=sort(d(d~=0)); %linearization of Walsh averages matrix
clear I J 
HLe=median(ld); %Hodges-Lehmann estimator
if n>15
    A=n*(n+1)/4; B=realsqrt(n*(n+1)*(2*n+1)/24); 
    Za=-realsqrt(2).*erfcinv(2.*0.975);
    T=fix(A-Za.*B);
else
    TC=[0 0 0 0 0 0 2 3 5 8 10 13 17 21 25];
    T=TC(n); clear TC
end
cintpem=ld([T+1 end-T]);
clear dff d ld T%clear unnecessary variable
fprintf('Median of differences (Hodges-Lehmann estimator): %0.4f\n',HLe)
fprintf('95%% Confidence interval: %0.4f  %0.4f\n',cintpem(1),cintpem(2))


%If the number of elements N<15 calculate the exact distribution of the
%signed ranks (the number of combinations is 2^N); else use the normal
%distribution approximation.

if n<=15
    ap=ff2n(n); %the all possibilities based on two-level full-factorial design.
    ap(ap~=1)=-1; %change 0 with -1
    k=1:1:n; 
    J=ap*k'; %all possible sums of ranks for k elements
    %to compute the p-value see how many values are more extreme of the observed
    %W and then divide for the total number of combinations
    p=length(J(abs(J)>=abs(W)))/length(J); %p-value
    %display results
    disp('The exact Wilcoxon distribution was used')
    disp(' ')
    fprintf('W\t\tp-value (2-tailed)\n')
    disp(tr)
    fprintf('%0.4f\t\t%0.4f\n',W,p)
    if nargout
        STATS.method='Exact distribution';
        STATS.W=W;
        STATS.p=p;
    end
else
    sW=sqrt((2*n^3+3*n^2+n-t)/6); %standard deviation
    zW=(abs(W)-0.5)/sW; %z-value with correction for continuity
    p=min([1 2*normcdf(zW)]); %p-value
    %display results
    disp('Sample size is good enough to use the normal distribution approximation')
    disp(' ')
    fprintf('W\t\tmW\t\tsW\t\tzW\t\tp-value (2-tailed)\n')
    disp(tr)
    fprintf('%0.4f\t0\t\t%0.4f\t\t%0.4f\t\t%0.15e\n',W,sW,zW,p)
    if nargout
        STATS.method='Normal approximation';
        STATS.W=W;
        STATS.mean=0;
        STATS.std=sW;
        STATS.z=zW;
        STATS.p=p;
    end
end
disp(tr)
disp(' ')

if plts
    xg=repmat([1 2],length(x1),1); yg=[x1; x2]';
    plot(xg,yg,'b.',xg',yg','r-')
    axis square
    set(gca,'XLim',[0 3],'XtickMode','manual','Xtick',0:1:3,'XtickLabel',{' ','Before','After',' '})
    title('Wilcoxon''s Plot')
    figure
    boxplot(yg(:),xg(:),'notch','on')
    set(gca,'XtickLabel',{'Before','After'})
end