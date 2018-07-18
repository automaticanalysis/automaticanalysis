function QuickCoReg(f1,f2,f3, opt)
%%% map f1 to f2

if nargin<4
    opt=1;
end

clear flags
flags.sep = [4 2]; % default
flags.params = [0 0 0  0 0 0];% default
flags.cost_fun = 'nmi'; % default; options are: 'mi'  - Mutual Information 'nmi' - Normalised Mutual Information 'ecc' - Entropy Correlation Coefficient  'ncc' - Normalised Cross Correlation
flags.tol = [0.02 0.02 0.02 0.001 0.001 0.001]; %default
flags.fwhm = [7 7]; % default
flags.graphics = ~spm('CmdLine') ;% default

VG = spm_vol(f2);
VF = spm_vol(f1);
x = spm_coreg(VG,VF,flags);
M  = spm_matrix(x);

save CoReg.mat x M;

%%% apply the estimated transformations to the header. look at spm_run_coreg_estwrite.m
M  = spm_matrix(x);
PO = {VF.fname};
MM = zeros(4,4,numel(PO));
MM = spm_get_space(VF.fname);
spm_get_space(VF.fname, M\MM);

%%% reslice the Mean CPS image to T1 space
clear P;
P{1}  = f2;
P{2}  = f1;
clear flags;
flags.mask   = 0;
flags.mean   = 0;
flags.interp = 1;
flags.which  = 1;
flags.wrap   = [0 0 0];
flags.prefix = 'r';

if opt==1
    spm_reslice(P,flags);
end


if nargin==3 && ~isempty(f3)
    [m h] = openIMG(f3);
    h2 = spm_vol(f1);
    h.mat = h2(1).mat;
    spm_write_vol(h,m);
    
    P{2} = f3;
    spm_reslice(P,flags);
end