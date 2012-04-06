% AA module - temporal filtering for EPIs
% Uses butterworth filter, and requires loading of entire series in one
%  go. This will strain the memory of some machines.
%  Rhodri Cusack Nov 2010

function [aap,resp]=aamod_highpassfilter_epi(aap,task,i,j)

resp='';

switch task
    case 'report'
        
    case 'doit'
        sesspath=aas_getsesspath(aap,i,j);
        imgs=aas_getimages_bystream(aap,i,j,'epi');
        nitem=0;
        nimg=size(imgs,1);
        for imgind=1:nimg
            V(imgind)=spm_vol(imgs(imgind,:));
            Y(imgind,:,:,:)=spm_read_vols(V(imgind));
        end;
        
        
        % Get TR from DICOM header
        if (length(aap.tasklist.currenttask.settings.TRs)==0)
            DICOMHEADERS=load(fullfile(sesspath,'dicom_headers'));
            aap.tasklist.currenttask.settings.TRs=DICOMHEADERS.DICOMHEADERS{1}.RepetitionTime/1000;
        end;

        % Regress out discrete cosine components to do filtering
        szY=size(Y);
        Y=reshape(Y,[nimg prod(szY(2:4))]);
        X0=spm_dctmtx(nimg,fix(2*(nimg*aap.tasklist.currenttask.settings.TRs)/aap.tasklist.currenttask.settings.HParam + 1));
        X0=X0(:,2:end);
        beta=X0\Y;
        Y=Y-X0*beta;
        Y=reshape(Y,szY);
        
        % Write the data
        aas_makedir(aap,sesspath);
        allimgs=[];
        for imgind=1:nimg
            [pth nme ext]=fileparts(V(imgind).fname);
            imgfn=fullfile(pth,['f' nme ext]);
            allimgs{imgind}=imgfn;
            V(imgind).fname=imgfn;
            spm_write_vol(V(imgind),squeeze(Y(imgind,:,:,:)));
        end;
        
        % Create the output stream
        aas_desc_outputs(aap,i,j,'epi',allimgs);
        
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;














