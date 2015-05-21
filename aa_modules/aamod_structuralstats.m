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
        streams = aas_getstreams(aap,'in');
        
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
        S.parts.desc = {'1=grey matter','2=white matter','3=csf'};
        S.TIV.mm3 = sum(S.parts.mm3);
        S.TIV.vox = sum(S.parts.vox);
        
        %% Corrected [RH]
        if aas_stream_has_contents(aap, subjind, 'seg8')
            jobs{1}.util{1}.tvol.matfiles = {aas_getfiles_bystream(aap,subjind,'seg8')};
            jobs{1}.util{1}.tvol.tmax = 3;
            jobs{1}.util{1}.tvol.mask = {fullfile(aap.directory_conventions.spmdir,'tpm','mask_ICV.nii,1')};
            jobs{1}.util{1}.tvol.outf = fullfile(aas_getsubjpath(aap,subjind),'structurals','SPM12_tissue_volumes');

            spm_jobman('initcfg')
            spm_jobman('run', jobs);
            
            tmp = textread(jobs{1}.util{1}.tvol.outf,'%s');
            tmp = textscan(tmp{2},'%s','Delimiter',',');
            for n=1:3
                S.parts.mm3_spm12(n) = str2num(tmp{1}{n+1})*1e6;
            end
            S.TIV.mm3_spm12 = sum(S.parts.mm3_spm12(1:3));
        end

        outfile = fullfile(fileparts(img), 'structuralstats.mat');
        save(outfile, 'S');
        aap = aas_desc_outputs(aap, subjind, 'structuralstats', outfile);
end