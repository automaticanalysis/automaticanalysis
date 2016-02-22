function STATS=testt(varargin)
%Student's t test for unpaired or paired samples.
% This file is applicable for equal or unequal sample sizes; for paired or 
% unpaired samples. When the test is unpaired, the Fisher-Snedecor F-test is
% performed to assess the equality of variance. If variances are not equal,
% Satterthwaite's approximate t test is performed.
% Testt requires powerStudent by Trujillo-Ortiz, A. and R. Hernandez-Walls. 
% URL http://www.mathworks.com/matlabcentral/fileexchange/2907
%
% Syntax: 	TESTT(X1,X2,TST,ALPHA,TAIL)
%      
%     Inputs:
%           X1 and X2 - data vectors (default = example data). 
%           TST - unpaired (0) or paired (1) test (default = 0).
%           ALPHA - significance level (default = 0.05).
%           TAIL - 1-tailed test (1) or 2-tailed test (2). (default = 2).
%     Outputs:
%           - t value.
%           - degrees of freedom.
%           - Confidence interval of means difference (for paired test)
%           - Critical value
%           - p-value
%
%      Example: 
%
%           X1=[77 79 79 80 80 81 81 81 81 82 82 82 82 83 83 84 84 84 84 85 ...
%           85 86 86 87 87];
% 
%           X2=[82 82 83 84 84 85 85 86 86 86 86 86 86 86 86 86 87 87 87 88 ...
%           88 88 89 90 90];
%
%           Calling on Matlab the function: testt
%
%           Answer is:
%
% FISHER-SNEDECOR F-TEST FOR EQUALITY OF VARIANCES
%  
% ------------------------------------------------------------
% F				DFn			DFd			p-value
% ------------------------------------------------------------
% 1.53788		24			24			0.29861
% ------------------------------------------------------------
% Variances are equal
% ------------------------------------------------------------
%  
% STUDENT'S T-TEST FOR UNPAIRED SAMPLES
%  
% ------------------------------------------------------------
% t				DF			  tail			p-value
% ------------------------------------------------------------
% 5.24110		48.0000			2			0.00000
% ------------------------------------------------------------
% It is a two-tailed hypothesis test.
% (The null hypothesis was statistically significative.)
%   
% Power is: 0.9989
%
% STATS=TESTT(...) returns a structure with all test(s) statistics
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2006). Student t-Test for unpaired or paired samples.
% http://www.mathworks.com/matlabcentral/fileexchange/12699

global n v alpha
%Input Error handling
args=cell(varargin);
nu=numel(args);
if isempty(nu) || nu==1
    error('Warning: Two data vectors are required')
elseif nu>5
    error('Warning: Max three input data are required')
end
default.values = {[77 79 79 80 80 81 81 81 81 82 82 82 82 83 83 84 84 84 84 85 ...
    85 86 86 87 87],[82 82 83 84 84 85 85 86 86 86 86 86 86 86 86 86 87 87 87 ...
    88 88 88 89 90 90],0,0.05,2};
default.values(1:nu) = args;
[x1 x2 tst alpha tail] = deal(default.values{:});
if ~isvector(x1) || ~isvector(x2)
   error('TESTT requires vector rather than matrix data.');
end 
if ~all(isfinite(x1)) || ~all(isnumeric(x1)) || ~all(isfinite(x2)) || ~all(isnumeric(x2))
    error('Warning: all X1 and X2 values must be numeric and finite')
end
if nu>2
    if ~isscalar(tst) || ~isfinite(tst) || ~isnumeric(tst) || isempty(tst)
        error('Warning: it is required a scalar, numeric and finite TST value.')
    end
    if tst ~= 0 && tst ~= 1 %check if tst is 0 or 1
        error('Warning: TST must be 0 for unpaired test or 1 for paired test.')
    end
    if tst==1
        if ((numel(x1) ~= numel(x2))),
            error('Warning: for paired test TESTT requires the data vectors to have the same number of elements.');
        end
    end
end
if nu>3
    if ~isscalar(alpha) || ~isnumeric(alpha) || ~isfinite(alpha) || isempty(alpha)
        error('Warning: it is required a numeric, finite and scalar ALPHA value.');
    end
    if alpha <= 0 || alpha >= 1 %check if alpha is between 0 and 1
        error('Warning: ALPHA must be comprised between 0 and 1.')
    end
end
if nu>4
    if ~isscalar(tail) || ~isfinite(tail) || ~isnumeric(tail) || isempty(tail)
        error('Warning: it is required a scalar, numeric and finite TAIL value.')
    end
    if tail ~= 2 && tail ~= 1 %check if tail is 1 or 2
        error('Warning: TAIL must be 1 or 2.')
    end
end
clear args default nu
tr=repmat('-',1,60);

switch tst
    case 0 %unpaired test
        n=[length(x1) length(x2)]; %samples sizes
        m=[mean(x1) mean(x2)]; %samples means
        v=[var(x1) var(x2)]; %samples variances
        %Fisher-Snedecor F-test
        if v(2)>v(1) 
            v=fliplr(v);
            m=fliplr(m);
            n=fliplr(n);
        end
        F=v(1)/v(2); %variances ratio
        DF=n-1;
        p = fcdf(F,DF(1),DF(2)); %p-value
        p = 2*min(p,1-p);
        if nargout
            STATS.Fvalue=F;
            STATS.DFn=DF(1);
            STATS.DFd=DF(2);
            STATS.FPvalue=p;
        end
        %display results
        disp('FISHER-SNEDECOR F-TEST FOR EQUALITY OF VARIANCES')
        disp(' ')
        disp(tr)
        fprintf('F\t\t\t\tDFn\t\t\tDFd\t\t\tp-value\n')
        disp(tr)
        fprintf('%0.5f\t\t\t%d\t\t\t%d\t\t\t%0.5f\n',F,DF,p)
        disp(tr)
        if p<alpha %unequal variances (Behrens-Welch problem)
            fprintf('Variances are different: Behrens-Welch problem\n')
            disp(tr)
            disp(' ')
            %Satterthwaite's approximate t test
            a=v./n; b=sum(a);
            denom=sqrt(b);
            gl=b^2/sum(a.^2./(n-1));
            disp('SATTERTHWAITE''S APPROXIMATE T-TEST FOR UNPAIRED SAMPLES')
            disp(' ')
            disp(tr)
        else %equal variances
            fprintf('Variances are equal\n')
            disp(tr)
            disp(' ')
            gl=sum(n)-2; %degrees of freedom
            s=sum((n-1).*v)/(sum(n)-2); %combined variance
            denom=sqrt(sum(s./n));
            disp('STUDENT''S T-TEST FOR UNPAIRED SAMPLES')
            disp(' ')
            disp(tr)
        end
        dm=diff(m); %Difference of means
        clear H n m v a b s %clear unnecessary variables
    case 1 %paired test
        disp('STUDENT''S T-TEST FOR PAIRED SAMPLES')
        disp(' ')
        disp(tr)
        n=length(x1); %samples size
        gl=n-1; %degrees of freedom
        d=x1-x2; %samples difference
        dm=mean(d); %mean of difference
        vc=tinv(1-alpha/tail,gl); %critical value
        ic=[abs(dm)-vc abs(dm)+vc]; %Confidence interval
        denom=sqrt((sum((d-dm).^2))/(n*(n-1))); %standard error of difference
        clear n d %clear unnecessary variables
        fprintf('Mean of difference\t\t\t\t')
        str=[num2str((1-alpha)*100) '%% C.I.\n'];
        fprintf(str)
        disp(tr)
        fprintf('%0.4f\t\t\t\t\t%0.4f\t\t\t%0.4f\n',abs(dm),ic)
        disp(tr)
end
t=abs(dm)/denom; %t value
p=(1-tcdf(t,gl))*tail; %t-value associated p-value
%display results
fprintf('t\t\t\t\tDF\t\t\t  tail\t\t\tp-value\n')
disp(tr)
fprintf('%0.5f\t\t\t%0.4f\t\t\t%d\t\t\t%0.5f\n',t,gl,tail,p)
disp(tr)
if nargout
    STATS.tvalue=t;
    STATS.tdf=gl;
    STATS.ttail=tail;
    STATS.tpvalue=p;
end
try
    powerStudent(t,gl,tail,alpha)
catch ME
    disp(ME)
    disp('I am trying to download the powerStudent function by Antonio Trujillo Ortiz from FEX')
    [F,Status]=urlwrite('http://www.mathworks.com/matlabcentral/fileexchange/2907-powerstudent?controller=file_infos&download=true','powerStudent.zip')
    if Status
        unzip(F)
        powerStudent(t,gl,tail,alpha)
    end
    clear F Status
end
