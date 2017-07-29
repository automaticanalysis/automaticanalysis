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
        bval = dlmread(bvals);
        bvec = dlmread(bvecs);
        [D11,D22,D33,D12,D13,D23,...
            W1111, W2222, W3333, W1112, W1113,...
            W1222, W2223, W1333, W2333, W1122,...
            W1133, W2233, W1123, W1223, W1233, dki.S0]=fun_DKI_ULLS_comp_rh(data_in,data_mask,bval,bvec);
        
        DT(:,:,:,6)=D23;
        DT(:,:,:,1)=D11;
        DT(:,:,:,2)=D22;
        DT(:,:,:,3)=D33;
        DT(:,:,:,4)=D12;
        DT(:,:,:,5)=D13;
        dki.DT = DT;
        
        KT(:,:,:,15)=W1233;
        KT(:,:,:,1)=W1111;
        KT(:,:,:,2)=W2222;
        KT(:,:,:,3)=W3333;
        KT(:,:,:,4)=W1112;
        KT(:,:,:,5)=W1113;
        KT(:,:,:,6)=W1222;
        KT(:,:,:,7)=W2223;
        KT(:,:,:,8)=W1333;
        KT(:,:,:,9)=W2333;
        KT(:,:,:,10)=W1122;
        KT(:,:,:,11)=W1133;
        KT(:,:,:,12)=W2233;
        KT(:,:,:,13)=W1123;
        KT(:,:,:,14)=W1223;
        dki.KT = KT;
        
        [dki.dMK, dki.dMD, dki.dS0]=fun_DKI_dMK_linear_rh(data_in,data_mask,bval);
        
        %% Calculate metrics
        [dki.MD, dki.FA, dki.AD, dki.RD, dki.L1, dki.L2, dki.L3, dki.V1, dki.V2, dki.V3] = fun_DTI_metrics(dki.DT,data_mask);
        [dki.MK, dki.AK, dki.RK, dki.AWF, dki.ADa, dki.ADe, dki.RDe, dki.tortu] = fun_DKI_metrics(dki.DT,dki.KT,data_mask);
                
        %% Now describe outputs
        V = spm_vol(betmask); V.dt = spm_type('float32');
        sesspath = aas_getpath_bydomain(aap,'diffusion_session',[subjind,diffsessind]);

        outstreams=aas_getstreams(aap,'output');        
        for outind=1:length(outstreams)
            metric = strrep(outstreams{outind},'dki_','');
            if ~isfield(dki,metric)
                aas_log(aap,false,sprintf('Metric %s for stream %s not exist!',metric,outstreams{outind}));
                continue; 
            end
            Y = dki.(metric);  
            nifti_write(fullfile(sesspath,[outstreams{outind} '.nii']),Y,outstreams{outind},V);
            aap=aas_desc_outputs(aap,'diffusion_session',[subjind,diffsessind],outstreams{outind},[outstreams{outind} '.nii']);
        end
        
end
end



