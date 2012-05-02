function [zthresh r2 r] = get_Z_thresh(N,NR,SigLevel)
% N = 20; %Sample Size (Number of Subjects for second level analysis, number of time points for first level analysis
% NR = 2; %Number of Regressors 2 for a Seed based analysis the seed plus a constant in the regressions equation.
% SigLevel = .001;
% N = length(SPM.xY.P);
% NR = size(SPM.xX.X,2);

 
df1 = NR-1;
df2 = N-df1-1;

Fs = 0:.001:100;
F_val = 1-cdf('F',Fs,df1,df2);
ind = find(F_val<SigLevel);
minF = Fs(ind(1));

 
r2 = 1/((1/(minF/df2*df1))+1);
r = sqrt(r2);
zthresh = atanh(r);
% [r2 r zthresh]

 