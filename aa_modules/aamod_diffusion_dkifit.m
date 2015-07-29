function [aap resp]=aamod_diffusion_dkifit(aap,task,subjind,diffsessind)
resp='';

switch task
    case 'report'
    case 'doit'
        %% Fetch inputs
        % Get nii filenames from stream
        diffinput=aas_getfiles_bystream(aap,'diffusion_session',[subjind diffsessind],'diffusion_data');
        bvals=aas_getfiles_bystream(aap,'diffusion_session',[subjind diffsessind],'bvals');
        bvecs=aas_getfiles_bystream(aap,'diffusion_session',[subjind diffsessind],'bvecs');
        betmask=aas_getfiles_bystream(aap,'diffusion_session',[subjind diffsessind],'BETmask');
        
        % Find which line of betmask contains the brain mask
        for betind=1:size(betmask,1)
            if strfind(betmask(betind,:),'bet_nodif_brain_mask')
                break
            end;
        end;        
        betmask=betmask(betind,:);
        
        %% Apply dkifit
        data_in = spm_read_vols(spm_vol(diffinput));
        data_mask = spm_read_vols(spm_vol(betmask));
        bval = importdata(bvals);
        bvec = importdata(bvecs);
        [dki.S0, dki.DT, dki.KT]=fun_DKI_ULLS_comp(data_in,data_mask,bval,bvec);
        [dki.dMK, dki.dMD, dki.dS0]=fun_DKI_dMK_linear(data_in,data_mask,bval);
        
        %% Calculate metrics
        [dki.MK, dki.AK, dki.RK] = fun_DKI_metrics(dki.DT,dki.KT,data_mask);
        [dki.MD, dki.FA, dki.AD, dki.RD, dki.L1, dki.L2, dki.L3, dki.V1, dki.V2, dki.V3] = fun_DTI_metrics(DT,data_mask);
        
        %% Now describe outputs
        V = spm_vol(betmask); V.dt = spm_type('float32');
        sesspath = aas_getpath_bydomain(aap,'diffusion_session',[subjind,diffsessind]);

        outstreams=aas_getstreams(aap,'out');        
        for outind=1:length(outstreams)
            metric = strrep(outstreams{outind},'dki_','');
            if ~exist(metric,'var')
                aas_log(aap,false,sprintf('Metric %s for stream %s not exist!',metric,outstreams{outind}));
                continue; 
            end
            Y = dki.(metric);  
            nifti_write(fullfile(sesspath,[outstreams{outind} '.nii']),Y,outstreams{outind},V);
            aap=aas_desc_outputs(aap,'diffusion_session',[subjind,diffsessind],outstreams{outind},[outstreams{outind} '.nii']);
        end
        
end
end



