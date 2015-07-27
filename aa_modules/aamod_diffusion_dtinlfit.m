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
        betmask=aas_getfiles_bystream(aap,'diffusion_session',[subjind diffsessind],'BETmask');
        
        % Find which line of betmask contains the brain mask
        for betind=1:size(betmask,1)
            if strfind(betmask(betind,:),'bet_nodif_brain_mask')
                break
            end;
        end;        
        betmask=betmask(betind,:);
        
        %% Apply dtinlfit
        data_in = spm_read_vols(spm_vol(diffinput));
        data_mask = spm_read_vols(spm_vol(betmask));
        bval = importdata(bvals);
        bvec = importdata(bvecs);
    
        if strcmp(aap.options.wheretoprocess, 'localsingle') && ~isempty(aap.directory_conventions.poolprofile) % local execution - use parfor for slices
            P = feval(aap.directory_conventions.poolprofile,size(data_mask,3));
            P.ResourceTemplate='-l nodes=^N^,mem=100MB,walltime=00:01:00';
            matlabpool(P);
            try
                parfor z = 1:size(data_in,3)
                    [S0(:,:,z), L1(:,:,z), L2(:,:,z), L3(:,:,z), V1(:,:,z,:), V2(:,:,z,:), V3(:,:,z,:)] = dti_slice(squeeze(data_in(:,:,z,:)), data_mask(:,:,z), bval, bvec);
                end
            catch ME
                fprintf('ERROR: %s\n', ME.message);
            end
            matlabpool close;
        else % cluster computing parfor is nor working
            for z = 1:size(data_in,3)
                [S0(:,:,z), L1(:,:,z), L2(:,:,z), L3(:,:,z), V1(:,:,z,:), V2(:,:,z,:), V3(:,:,z,:)] = dti_slice(squeeze(data_in(:,:,z,:)), data_mask(:,:,z), bval, bvec);
            end
        end
        
        AD = L1;
        RD = (L2+L3)/2;
        MD = (L1+L2+L3)/3;
        FA = sqrt(3/2)*...
            sqrt(var(cat(4, L1, L2, L3),[],4)*2./...
            (L1.^2+L2.^2+L3.^2));
        
        % Now describe outputs
        V = spm_vol(betmask); V.dt = spm_type('float32');
        sesspath = aas_getpath_bydomain(aap,'diffusion_session',[subjind,diffsessind]);

        outstreams=aap.tasklist.currenttask.outputstreams;        
        for outind=1:length(outstreams.stream)
            Y = eval(strrep(outstreams.stream{outind},'dti_',''));
            Y(isnan(Y)) = 0;
            if ndims(Y) == 3
                nifti_write(fullfile(sesspath,[outstreams.stream{outind} '.nii']),...
                    Y, outstreams.stream{outind},V);
            else % 4D
                for z = 1:size(Y,4)
                    outfiles{z} = fullfile(sesspath,sprintf('%s-%04d.nii',outstreams.stream{outind},z));
                    nifti_write(outfiles{z},Y(:,:,:,z),outstreams.stream{outind},V);
                end
                spm_file_merge(char(outfiles),fullfile(sesspath,[outstreams.stream{outind} '.nii']),0);
                delete(fullfile(sesspath,[outstreams.stream{outind} '-*']));
            end
            aap=aas_desc_outputs(aap,'diffusion_session',[subjind,diffsessind],outstreams.stream{outind},[outstreams.stream{outind} '.nii']);
        end;
        
end
end

%% DTI wrapper

function [S0, L1, L2, L3, V1, V2, V3] = dti_slice(data, mask, b, bv)

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
        
        [S0(x,y), DT] = dti_nlfit(vdata, b, bv);
        
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