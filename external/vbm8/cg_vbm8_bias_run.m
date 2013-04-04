function out = cg_vbm8_bias_run(job)
% Correct bias between an image pair
%__________________________________________________________________________
% Christian Gaser
% $Id: cg_vbm8_bias_run.m 405 2011-04-12 11:31:21Z gaser $

for i=1:numel(job.subj),
    out(i).files = cell(numel(job.subj(i).mov)-1,1);
    for j=2:numel(job.subj(i).mov)
        [pth,nam,ext,num] = spm_fileparts(job.subj(i).mov{j});
        out(i).files{j-1} = fullfile(pth,['m', nam, ext, num]);
        run_bias_correction(job.subj(i).mov{j},job.subj(i).mov{1},job.bias_opts);
    end
end;
%_______________________________________________________________________

%_______________________________________________________________________
function run_bias_correction(PF,PG,bias_opts)

VG  = spm_vol(PG);
VF  = spm_vol(PF);

bo  = bias_opts;
if isfinite(bo.fwhm) & (bo.fwhm > 0)
    dat = bias_correction(VF,VG,bo.nits,bo.fwhm,bo.reg,bo.lmreg);
else
    dat = spm_read_vols(VF);
end

[pth,nam,ext,num] = spm_fileparts(VF.fname);
VF.fname = fullfile(pth,['m', nam, ext, num]);
VF.descrip = 'Bias corrected image';

spm_write_vol(VF, dat);

return;
%_______________________________________________________________________

%_______________________________________________________________________
function t = transf(B1,B2,B3,T)
d2 = [size(T) 1];
t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
t  = B1*t1*B2';
return;
%_______________________________________________________________________

%_______________________________________________________________________
function dat = bias_correction(VF,VG,nits,fwhm,reg1,reg2)
% This function is intended for doing bias correction prior
% to longitudinal VBM.  
% A version of the first image is returned out, that has the
% same bias as that of the second image.  

vx      = sqrt(sum(VF.mat(1:3,1:3).^2));
d       = VF.dim(1:3);
sd      = vx(1)*VF.dim(1)/fwhm; d3(1) = ceil(sd*2); krn_x   = exp(-(0:(d3(1)-1)).^2/sd.^2)/sqrt(vx(1));
sd      = vx(2)*VF.dim(2)/fwhm; d3(2) = ceil(sd*2); krn_y   = exp(-(0:(d3(2)-1)).^2/sd.^2)/sqrt(vx(2));
sd      = vx(3)*VF.dim(3)/fwhm; d3(3) = ceil(sd*2); krn_z   = exp(-(0:(d3(3)-1)).^2/sd.^2)/sqrt(vx(3));
Cbias   = kron(krn_z,kron(krn_y,krn_x)).^(-2)*prod(d)*reg1;
Cbias   = sparse(1:length(Cbias),1:length(Cbias),Cbias,length(Cbias),length(Cbias));
B3bias  = spm_dctmtx(d(3),d3(3));
B2bias  = spm_dctmtx(d(2),d3(2));
B1bias  = spm_dctmtx(d(1),d3(1));
lmRb    = speye(size(Cbias))*prod(d)*reg2;
Tbias   = zeros(d3);

% correct for different means
globalF = spm_global(VF);
globalG = spm_global(VG);
VF.pinfo(1:2,:) = globalG/globalF*VF.pinfo(1:2);

ll = Inf;
try
    spm_plot_convergence('Init','Bias Correction','- Log-likelihood','Iteration');
catch
    spm_chi2_plot('Init','Bias Correction','- Log-likelihood','Iteration');
end
for subit=1:nits,

    % Compute objective function and its 1st and second derivatives
    Alpha = zeros(prod(d3),prod(d3)); % Second derivatives
    Beta  = zeros(prod(d3),1); % First derivatives
    oll   = ll;
    ll    = 0.5*Tbias(:)'*Cbias*Tbias(:);

    for z=1:VF.dim(3),
        M1  = spm_matrix([0 0 z]);
        M2  = VG.mat\VF.mat*M1;
        f1o = spm_slice_vol(VF,M1,VF.dim(1:2),0);
        f2o = spm_slice_vol(VG,M2,VF.dim(1:2),0);
        f1o(~isfinite(f1o)) = 0;
        f2o(~isfinite(f2o)) = 0;
        msk = (f1o==0) & (f2o==0);
        f1o(msk) = 0;
        f2o(msk) = 0;
        ro       = transf(B1bias,B2bias,B3bias(z,:),Tbias);
        msk      = abs(ro)>0.01; % false(d(1:2));

        % Use the form based on an integral for bias that is
        % far from uniform.
        f1  = f1o(msk);
        f2  = f2o(msk);
        r   = ro(msk);
        e   = exp(r);
        t1  = (f2.*e-f1);
        t2  = (f1./e-f2);
        ll  = ll + 1/4*sum(sum((t1.^2-t2.^2)./r));
        wt1 = zeros(size(f1o));
        wt2 = zeros(size(f1o));
        wt1(msk) = (2*(t1.*f2.*e+t2.*f1./e)./r + (t2.^2-t1.^2)./r.^2)/4;
        wt2(msk) = ((f2.^2.*e.^2-f1.^2./e.^2+t1.*f2.*e-t2.*f1./e)./r/2 ...
                 - (t1.*f2.*e+t2.*f1./e)./r.^2 + (t1.^2-t2.^2)./r.^3/2);

        % Use the simple symmetric form for bias close to uniform
        f1  = f1o(~msk);
        f2  = f2o(~msk);
        r   = ro(~msk);
        e   = exp(r);
        t1  = (f2.*e-f1);
        t2  = (f1./e-f2);
        ll  = ll + (sum(t1.^2)+sum(t2.^2))/4;
        wt1(~msk) = (t1.*f2.*e-t2.*f1./e)/2;
        wt2(~msk) = ((f2.*e).^2+t1.*f2.*e + (f1./e).^2+t2.*f1./e)/2;

        b3    = B3bias(z,:)';
        Beta  = Beta  + kron(b3,spm_krutil(wt1,B1bias,B2bias,0));
        Alpha = Alpha + kron(b3*b3',spm_krutil(wt2,B1bias,B2bias,1));
    end;
    try
        spm_plot_convergence('Set',ll/prod(d));
    catch
        spm_chi2_plot('Set',ll/prod(d));
    end

    if subit > 1 && ll>oll,
        % Hasn't improved, so go back to previous solution
        Tbias = oTbias;
        ll    = oll;
        lmRb  = lmRb*10;
    else
        % Accept new solution
        oTbias = Tbias;
        Tbias  = Tbias(:);
        Tbias  = Tbias - (Alpha + Cbias + lmRb)\(Beta + Cbias*Tbias);
        Tbias  = reshape(Tbias,d3);
    end;
end;

dat = zeros(VF.dim(1:3));

for z=1:VF.dim(3),
    M1 = spm_matrix([0 0 z]);
    tmp = spm_slice_vol(VF,M1,VF.dim(1:2),0);
    tmp(~isfinite(tmp)) = 0;
    dat(:,:,z) = tmp;
    r  = transf(B1bias,B2bias,B3bias(z,:),Tbias);
    dat(:,:,z) = dat(:,:,z)./exp(r);
end;

return;
%_______________________________________________________________________
