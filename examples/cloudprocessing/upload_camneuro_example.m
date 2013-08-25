if (~exist('cc','var') || isempty(cc))
    aa_ver4_devel
    import camneuro.aacc.*;
    cc=Aacc();
    cc.login();
end;

% CHOOSE ONE OR MORE OF
% (1) SPECIFY SINGLE SUBJECTS
% Replace with paths to raw data directories...
%cc.syncDicom('/mridata/cbu/CBU0123456_MR012345');
% Add as many lines as you like...
%cc.syncDicom('/mridata/cbu/CBU0123457_MR012345');
%cc.syncDicom('/mridata/cbu/CBU0123458_MR012345');

% (2) UPLOAD EVERYTHING WITH A GIVEN MRxxxxx NUMBER

pth='/mridata/cbu';
fn=dir(fullfile(pth,'*_MR10023'));
for fnind=1:length(fn)
    cc.syncDicom(fullfile(pth,fn(fnind).name));
end;
