% AA module - temporal filtering for EPIs
% Uses SPM style filter, and requires loading of entire series in one
%  go. This will strain the memory of some machines.
%  Rhodri Cusack Nov 2010

function [aap,resp]=aamod_highpassfilter_epi(aap,task,i,j)

resp='';

switch task
    case 'report'
        
    case 'doit'
        sesspath=aas_getsesspath(aap,i,j);
        imgs=aas_getfiles_bystream(aap,i,j,'epi');
        nitem=0;
        
        % Read in data, may be 4D or compressed 
        [V, Y , ~, nimg]=aas_spm_vol(imgs);
        
        
        % Get TR from DICOM header
        if isempty(aap.tasklist.currenttask.settings.TRs)
            DICOMHEADERS=load(aas_getfiles_bystream(aap,i,j,'epi_dicom_header'));
            aap.tasklist.currenttask.settings.TRs=DICOMHEADERS.DICOMHEADERS{1}.RepetitionTime/1000;
        end;
        
        % Regress out discrete cosine components to do filtering
        
        Y=permute(Y,[4 1 2 3]);  % Time becomes first dimension, temporarily
        szY=size(Y);
        Y=reshape(Y,[nimg prod(szY(2:4))]);
        X0=spm_dctmtx(nimg,fix(2*(nimg*aap.tasklist.currenttask.settings.TRs)/aap.tasklist.currenttask.settings.HParam + 1));
        X0=X0(:,2:end);
        beta=X0\Y;
        Y=Y-X0*beta;
        Y=reshape(Y,szY);
     
        Y=permute(Y,[2 3 4 1]);   % Put time back at end
        
        % Write the data, now as 4D
        aas_makedir(aap,sesspath);
        allimgs=[];
        
        outfn=V(1).fname;
        [pth nme ext]=fileparts(outfn);
        if ~strcmp(nme(end-2:end),'_4D')
            nme=[nme '_4D'];
        end;
        outfn=fullfile(pth,['f' nme ext]);
        
        
        % Write file out a single 4D nii
        aas_spm_write_vol(V(1),Y,outfn);
        
        % Create the output stream
        aap = aas_desc_outputs(aap,i,j,'epi',outfn);
        
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
