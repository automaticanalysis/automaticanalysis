%

function [aap resp]=aamod_probtrackxsummarize_group(aap,task)

resp='';

switch task
    case 'report'
    case 'doit'
        fslext=aas_getfslext(aap);
        sessind=1;  
        % Load seed ROI in individual diffusion space
        maskfn=aas_getfiles_bystream(aap,'study',[],'tractography_seeds_mni');
        % Reslice into space of seeds
        s2t=aas_getfiles_bystream(aap,'diffusion_session',[1 sessind],'seeds_to_mni_space');
        spm_reslice(strvcat(s2t,maskfn));
        [pth nme ext]=aas_fileparts(maskfn);
        maskfn_r=fullfile(pth,['r' nme ext]);
        
        [Vmask Ymask allXYZ]=aas_spm_vol(maskfn_r);
        Ymask=Ymask(:)>0.1;
        maskind=find(Ymask(:));
        seeddim=sum(Ymask(:));
        
              
        rowoffset=0;
        nsubj=length(aap.acq_details.subjects);
        for subjind=1:nsubj
            s2t=aas_getfiles_bystream(aap,'diffusion_session',[subjind sessind],'seeds_to_mni_space');
            if subjind==1
                ntarg=size(s2t,1);
                seed2target=zeros(ntarg*nsubj,seeddim);
            else
                if size(s2t,1)~=ntarg
                    aas_log(aap,true,sprintf('All seeds_to_mni_space must have same number of targets but subject %d differs.',subjind));
                end;
            end;
            
            for targind=1:ntarg
                [V Y]=aas_spm_vol(s2t(targind,:));
                seed2target(targind+rowoffset,:)=Y(Ymask);
            end;
            rowoffset=rowoffset+ntarg;
        end;
        
        s2tfn=fullfile(aas_getstudypath(aap),'group_tractography_seed2target.mat');
        save(s2tfn,'seed2target','maskind','allXYZ','maskfn');
        aap=aas_desc_outputs(aap,'study',[],'seed2target',s2tfn);
end

