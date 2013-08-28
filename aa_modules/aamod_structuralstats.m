function [aap,resp] = aamod_structuralstats(aap, task, subjind)
% AAMOD_STRUCTURALSTATS Get total tissue class volumes and TIV for segmented images.
%
% Get volume in mm^3 of GM, WM, CSF, and total intracranial volume (TIV).
%
% input streams:    native_grey
%                   native_white
%                   native_csf
%
% output stream:    structuralstats

resp='';

switch task
    case 'domain'
        resp = 'subject';
    case 'report'
        resp = 'Get statistics on segmented structural images.'
    case 'doit'

        streams = {'native_grey', 'native_white', 'native_csf'};

        for k=1:3
            img = aas_getfiles_bystream(aap, subjind, streams{k});

            V = spm_vol(img);
            Y = spm_read_vols(V);
            spacedesc = spm_imatrix(V.mat);
            volume = prod(abs(spacedesc(7:9)));
            Y = Y(~isnan(Y));
            S.parts.vox(k) = sum(Y(:));
            S.parts.mm3(k) = S.parts.vox(k)*volume;
        end % going through tissue classes

        S.parts.desc = {'1=grey matter','2=white matter','3=csf'};
        S.TIV.mm3 = sum(S.parts.mm3);
        S.TIV.vox = sum(S.parts.vox);

        [pth, nm, ext] = fileparts(img);

        outfile = fullfile(pth, 'structuralstats.mat');
        save(outfile, 'S');
        aap = aas_desc_outputs(aap, subjind, 'structuralstats', outfile);
end