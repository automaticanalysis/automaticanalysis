% aa module - calculates spectrum in time of different tissue compartments
% Rhodri Cusack BMI Western Nov 2013

function [aap,resp]=aamod_tissue_wavelets(aap,task,subj,sess)
resp='';

switch task
    case 'report'
        
    case 'doit'
        
        tissues= aap.tasklist.currenttask.settings.tissuelist.tissue;
        dwtlevels=aap.tasklist.currenttask.settings.dwtlevels;
        dwtname=aap.tasklist.currenttask.settings.dwtname; 

        % EPIs
        epi_fns = aas_getfiles_bystream(aap,'session',[subj,sess],'epi');
        
        % And tissue segementations
        seg_fns={};
        for tissind=1:length(tissues)
            % Get filename of this compartment
            seg_fns{tissind}= aas_getfiles_bystream(aap,'subject',subj,['tissue_' tissues{tissind}]);
            
            % Reslice to space of EPIs
            spm_reslice(strvcat(epi_fns(1,:),seg_fns{tissind}),struct('which',1));
            
            % Load up resliced image
            [pth nme ext]=aas_fileparts(seg_fns{tissind});
            seg_V{tissind}=spm_vol(fullfile(pth,['r' nme ext]));
            mask=spm_read_vols(seg_V{tissind})>0.5;
            tissmask(:,tissind)=mask(:);
        end;
        
        % Load up EPIs
        
        tiss_timecourse=zeros([size(tissmask,4) size(epi_fns,1)]);

        % 4D?
        if size(epi_fns,1)==1
            Vepi=spm_vol(epi_fns);
            Yepi_timecourse=spm_read_vols(Vepi);
        else
            % no 3D
            for fnind=1:size(epi_fns,1)
                Vepi=spm_vol(epi_fns(fnind,:));

                Yepi=spm_read_vols(Vepi);
                Yepi_timecourse(:,fnind)=Yepi(:);
            end;
        end;
        
         % For first analysis, lets collapse across voxels within each of
        % the three tissue compartments
        dim=size(Yepi_timecourse);
        
        % Reshape 
        Yepi_timecourse=reshape(Yepi_timecourse,[prod(dim(1:3)) dim(4)]);
        % Rescale so mean=100
        voxmean=mean(Yepi_timecourse,2);
        mask=voxmean~=0;
        Yepi_timecourse(mask,:)=100*Yepi_timecourse(mask,:)./repmat(voxmean(mask),[1 dim(4)]);
        
        sesspth=aas_getsesspath(aap,subj,sess);
        tiss_timecourse_fn=fullfile(sesspth,'tiss_timecourse.mat');
        save(tiss_timecourse_fn,'tiss_timecourse');
        aap=aas_desc_outputs(aap,'session',[subj sess],'tiss_timecourse',tiss_timecourse_fn);
        
        % Discrete wavelet transform
        tiss_wavelet=cell(length(tissues),1);
        for tissind=1:length(tissues)
            dat=Yepi_timecourse(tissmask(:,tissind),:);
            tiss_wavelet{tissind}=mdwtdec('r',dat,dwtlevels,dwtname);   
        end;

        tiss_powerspect_fn=fullfile(sesspth,'tiss_wavelets.mat');
        save(tiss_powerspect_fn,'tiss_wavelet','tissmask','tissues');
        aap=aas_desc_outputs(aap,'session',[subj sess],'tiss_wavelets',tiss_powerspect_fn);
  
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;



