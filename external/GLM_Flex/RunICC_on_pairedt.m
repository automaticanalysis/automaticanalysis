function RunICC_on_pairedt
% keyboard;
load oo
load I;

dfJ = 2;
dfT = numel(I.Scans)/2;


F = I.F;
Y = oo.Y{1};

X = I.X;
X(isnan(X))=0;

beta = oo.beta1{1};

[r c cc rr] = MakeContrastMatrix('a',F.XX(:,3:end),[dfT],F.XX,X);
BSS = LoopEstimate(beta,X,r);
[r c cc rr] = MakeContrastMatrix('a',F.XX(:,1:2),[dfJ],F.XX,X);
WSS = LoopEstimate(beta,X,r);

er = eye(size(X,1))-(X*pinv(X));
RE = LoopEstimate(Y,1,er);
[z1 z2 z3] = svd(X);
tol = max(size(X))*max(abs(diag(z2)))*eps;
df1 = sum(diag(z2)>tol);
df1 = size(X,1)-df1;
df = df1;

MSE = RE/df;

dat = ((BSS/(dfT-1)) - MSE) ./ ((BSS/(dfT-1)) + MSE + (dfJ*((WSS/(dfJ-1))-MSE)/dfT));

[m,h] = openIMG('NN.nii');
h.fname = 'ICC.nii';
m(oo.vec{1}) = dat;
spm_write_vol(h,m);

%%
% % % dat = ((BSS/(dfT)) - (WSS/(dfJ))) ./ ((BSS/(dfT)) - ((dfT-1)*(WSS/(dfJ))));
% % % - ((dfT-1)*(WSS/(dfJ-1))));
% % 
% % % tmp = openIMG('ResMS_01.nii');
% dat = ((BSS/99)-(WSS/1))./((BSS/99) +(WSS/1));
% % dat = ((WSS/1)-(BSS/99))./((WSS/1)+(BSS));
% % 
% % % figure(20); clf; hist(dat,100); shg
% % % 
% % % %%
% % % dat = ((BSS/99)-(WSS/1))
% [m,h] = openIMG('NN.nii');
% h.fname = 'ICC1.nii';
% m(oo.vec{1}) = dat;
% spm_write_vol(h,m);


