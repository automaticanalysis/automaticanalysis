function [aap,resp]=aamod_normalisebytotalgrey(aap, task, subjind)
%AAMOD_NORMALISEBYTOTALGREY Divide voxelwise GM by TGM as proportional scaling for VBM.
%
% Divide structural images by total gray matter volume as a form of
% proportional scaling. Output images have a 'g' prepended (for global scaling).
%
% Requires aamod_structuralstats to be run first to get total gray
% matter.
%
% Input streams:    structuralstats
%                   normalised_grey
%
% Output stream:    tgmscaled_normalised_grey
%
% For more about global scaling in VBM, see:
%
%     Peelle JE, Cusack R, Henson RNA (2012) Adjusting for global effects in voxel-based
%     morphometry: Gray matter decline in normal aging. NeuroImage 60:1503-1516.

resp='';

% possible tasks 'doit','report','checkrequirements'

switch task
    case 'domain'
        resp='subject';
    case 'report'
        resp='Normalise images by total grey matter.';

    case 'doit'

        % get files
        structstats = aas_getfiles_bystream(aap, subjind, 'structuralstats');
        gmimg = aas_getfiles_bystream(aap, subjind,  'normalised_grey');

        % get the total grey matter (TGM)
        load(structstats);
        tgm = S.parts.mm3(1) * 1000; % to get litres instead of mL - makes  display a little easier

        % perform division and write out new image
        V = spm_vol(gmimg);
        Y = spm_read_vols(V);
        Y = Y./tgm;
        [pth, nm, ext] = fileparts(V.fname);
        V.fname = fullfile(pth, ['g' nm ext]);
        spm_write_vol(V, Y);

        % describe outputs
        aap = aas_desc_outputs(aap, subjind, 'tgmscaled_normalised_grey', V.fname);
end