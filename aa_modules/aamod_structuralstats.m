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
        resp = 'Get statistics on segmented structural images.';
    case 'doit'
        streams = aas_getstreams(aap,'input');
        S.parts.desc = {'1=grey matter','2=white matter','3=csf'};
        
        %% Conventional
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
        S.TIV.mm3 = sum(S.parts.mm3);
        S.TIV.vox = sum(S.parts.vox);
        
        %% Corrected [RH]
        if aas_stream_has_contents(aap, subjind, 'seg8')
            job.matfiles = {aas_getfiles_bystream(aap,subjind,'seg8')};
            job.tmax = 3;
            job.mask = {fullfile(aap.directory_conventions.spmdir,'tpm','mask_ICV.nii,1')};
            job.outf = '';

            out = spm_run_tissue_volumes('exec', job);
            
            S.parts.mm3_spm12 = [out.vol1 out.vol2 out.vol3]*1e6;
            S.TIV.mm3_spm12 = sum(S.parts.mm3_spm12(1:3));
        end

        outfile = fullfile(fileparts(img), 'structuralstats.mat');
        save(outfile, 'S');
        aap = aas_desc_outputs(aap, subjind, 'structuralstats', outfile);
end