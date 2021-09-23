% AA module: re-write EPIs with float precision

function [aap,resp]=aamod_make_epis_float(aap,task,subj,sess)

resp='';

switch task
    case 'doit'
        imfn = aas_getimages_bystream(aap,subj,sess,'epi');
        V = spm_vol(imfn);
        nvol = length(V);
        xyz = spm_read_vols(V);
        for v = 1:nvol
          V(v).fname = imfn(v,:);
          V(v).dt = [spm_type('float64') spm_platform('bigend')];
          spm_write_vol(V(v),xyz(:,:,:,v));
        end
        aap=aas_desc_outputs(aap,subj,sess,'epi',imfn);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
