function [aap resp]=aamod_diffusion_dtinlfit(aap,task,subjind,diffsessind)
resp='';

switch task
    case 'report'
    case 'doit'
        %% Fetch inputs
        % Get nii filenames from stream
        diffinput=aas_getfiles_bystream(aap,'diffusion_session',[subjind diffsessind],'diffusion_data');
        bvals=aas_getfiles_bystream(aap,'diffusion_session',[subjind diffsessind],'bvals');
        bvecs=aas_getfiles_bystream(aap,'diffusion_session',[subjind diffsessind],'bvecs');
        betmask=cellstr(aas_getfiles_bystream(aap,'diffusion_session',[subjind diffsessind],'BETmask'));
        
        % Find which line of betmask contains the brain mask
        betmask=betmask{cellfun(@(x) ~isempty(regexp(x,'bet_.*nodif_brain_mask', 'once')), betmask)};
        
        %% Apply dtinlfit
        data_in = spm_read_vols(spm_vol(diffinput));
        data_mask = spm_read_vols(spm_vol(betmask));
        bval = dlmread(bvals);
        bvec = dlmread(bvecs);
    
        [D11,D22,D33,D12,D13,D23, dti.S0]=fun_DTI_UNLS_comp_rh(data_in,data_mask,bval,bvec);
        DT(:,:,:,6)=D23;
        DT(:,:,:,1)=D11;
        DT(:,:,:,2)=D22;
        DT(:,:,:,3)=D33;
        DT(:,:,:,4)=D12;
        DT(:,:,:,5)=D13;
        [dti.MD, dti.FA, dti.AD, dti.RD, dti.L1, dti.L2, dti.L3, dti.V1, dti.V2, dti.V3] = fun_DTI_metrics(DT,data_mask);
              
        % Now describe outputs
        V = spm_vol(betmask); V.dt = spm_type('float32');
        sesspath = aas_getpath_bydomain(aap,'diffusion_session',[subjind,diffsessind]);

        outstreams=aas_getstreams(aap,'output');        
        for outind=1:length(outstreams)
            metric = strrep(outstreams{outind},'dti_','');
            if ~isfield(dti,metric)
                aas_log(aap,false,sprintf('Metric %s for stream %s not exist!',metric,outstreams{outind}));
                continue; 
            end
            Y = dti.(metric);  
            nifti_write(fullfile(sesspath,[outstreams{outind} '.nii']),Y,outstreams{outind},V);
            aap=aas_desc_outputs(aap,'diffusion_session',[subjind,diffsessind],outstreams{outind},[outstreams{outind} '.nii']);
        end
end
end

%% DTI wrapper

function [S0, L1, L2, L3, V1, V2, V3] = dti_slice(data, mask, b, bv, tol_b0)

S0 = mask*0;
L1 = mask*0;
L2 = mask*0;
L3 = mask*0;
V1 = repmat(mask*0,[1 1 3]);
V2 = repmat(mask*0,[1 1 3]);
V3 = repmat(mask*0,[1 1 3]);

for x = 1:size(data,1)
    for y = 1:size(data,2)
        
        vdata = squeeze(data(x,y,:))';
        
        if ~mask(x,y), continue; end
        
        [S0(x,y), DT] = dti_nlfit(vdata, b, bv, tol_b0);
        
        [V,D] = eig(DT);

        L1(x,y) = D(3,3);
        L2(x,y) = D(2,2);
        L3(x,y) = D(1,1);
        
        V1(x,y,:) = V(:,3);
        V2(x,y,:) = V(:,2);
        V3(x,y,:) = V(:,1);
    end
end
end