function [xBF] = aas_get_custom_BF(aap, xBF, TR)
% Use this function to create or use custom basis functions in SPM. We
% will assume that the custom BFs are stored in the aas_custom_BFs.mat data
% structure. Create a new case in the switch statement here to do what you
% need to. 

dt = TR/xBF.T;

if ~exist('aa_custom_BFs.mat'), aas_log(aap, 1, 'Can''t find custom BF file!'); end
load('aa_custom_BFs.mat');

bfI = find(strcmp(xBF.name, {aa_custom_BFs.name}));
if isempty(bfI), aas_log(aap, 1, sprintf('Invalid basis function name: %s', xBF.name)); end

hrfTR = aa_custom_BFs(bfI).fs;
hrfLen = size(aa_custom_BFs(bfI).xBF, 1) * hrfTR;

Thrf = 0 : hrfTR : hrfLen;
Thrf = Thrf(1:end-1);

Tdes = 0 : dt : hrfLen;
Tdes = Tdes(1:end-1);

bf = [];
for b = 1 : size(aa_custom_BFs(bfI).xBF, 2)
    bf = [bf interp1(Thrf, aa_custom_BFs(bfI).xBF(:,b), Tdes, 'linear', 'extrap')'];
end


% Orthogonalise and fill in basis function structure
%--------------------------------------------------------------------------
xBF.bf  =  spm_orth(bf);
