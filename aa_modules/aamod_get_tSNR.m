%AA module - calculate tSNR and get distribution 

function [aap,resp]=aamod_get_tSNR(aap,task,subj,sess)

resp = '';
 
switch task
    case 'report'
        
    case 'doit'
        
        imgs=aas_getimages_bystream(aap,subj,sess,'epi');
        
        sesspath=aas_getsesspath(aap,subj,sess);
 
        Yepi=[];
        % Load up all volumes
        for ind=1:size(imgs,1)
            Vepi=spm_vol(sprintf('%s',imgs(ind,:)));
            Yepi(:,:,:,ind)=spm_read_vols(Vepi);
        end;
        
        % calculate mean and std through time
        Ymn(:,:,:)=mean(Yepi,4);
        Ystd(:,:,:)=std(Yepi,[],4);
        
        % Load TR
        dcmfn=aas_getimages_bystream(aap,subj,sess,'epi_dicom_header');
        H=load(dcmfn);
        TR=H.DICOMHEADERS{1}.RepetitionTime;
        TE=H.DICOMHEADERS{1}.EchoTime; 
        
        % Calculate tSNR
        Ysnr=TE*Ymn./Ystd/sqrt(TR);
        
        % write out SNR map 
        Vsnr=Vepi; 
        Vsnr.fname=fullfile(sesspath,'tSNR_scaled.nii');
        spm_write_vol(Vsnr,Ysnr);
        
        %Find mask
        thresh=aap.tasklist.currenttask.settings.mask_threshold;
        Ymeanepi_thresh=Ymn>thresh;
        
        %write out mean EPI
        Vmn=Vepi; 
        Vmn.fname=fullfile(sesspath,'meanEPI.nii');
        spm_write_vol(Vmn,Ymn);
        
        %write out mask 
        Vmask=Vepi; Ymask=Ymeanepi_thresh; 
        Vmask.fname=fullfile(sesspath,'meanEPI_mask.nii');
        spm_write_vol(Vmask,Ymask);
        
        %get histogram of SNR within mask 
        dat=Ysnr(:,:,:);        
        dat=dat(Ymeanepi_thresh(:)>0);
        dat=dat(~isnan(dat) & ~isinf(dat));
        [N, X]=hist(dat,[1:4:max(dat)+50]);
        [maxN, indN]=max(N);
        
        fprintf('mode of tSNR is %d\n',X(indN)); 
        
        h=figure(12);
        plot(X,N);
        xlabel('tSNR');
        ylabel('number of voxels');
        fprintf('\n');
             
        save(fullfile(sesspath,'tSNR_hist.mat'),'dat','N','X','maxN','indN');  
        saveas(h,fullfile(sesspath,'tSNR_dist'),'png'); 
end

