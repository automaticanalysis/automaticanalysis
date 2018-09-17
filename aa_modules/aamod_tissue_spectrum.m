% aa module - calculates spectrum in time of different tissue compartments
% Rhodri Cusack BMI Western Nov 2013

function [aap,resp]=aamod_tissue_spectrum(aap,task,subj,sess)
resp='';

switch task
    case 'report'
        
    case 'doit'
        
        tissues= aap.tasklist.currenttask.settings.tissuelist.tissue;
        if ~iscell(tissues)
            tissues={tissues};
        end;
        
        % EPIs
        epi_fns = aas_getfiles_bystream(aap,'session',[subj,sess],'epi');
        
        % dicom_header
        TR=double(aap.tasklist.currenttask.settings.TR);
        if isempty(TR)
            epi_header= aas_getfiles_bystream(aap,'session',[subj,sess],'epi_header');
            load(epi_header,'DICOMHEADERS');
            TR=DICOMHEADERS{1}.RepetitionTime/1000;
        end;
        
        % And tissue segementations
        seg_fns={};
        tiss_mask=[];
        tiss_weight=[];
        for tissind=1:length(tissues)
            % Get filename of this compartment
            seg_fns{tissind}= aas_getfiles_bystream(aap,'subject',subj,['tissue_' tissues{tissind}]);
            
            % Reslice to space of EPIs
            spm_reslice(strvcat(epi_fns(1,:),seg_fns{tissind}),struct('which',1));
            
            % Load up resliced image
            [pth nme ext]=aas_fileparts(seg_fns{tissind});
            seg_V{tissind}=spm_vol(fullfile(pth,['r' nme ext]));
            mask=spm_read_vols(seg_V{tissind})>0.5;
            tiss_mask(:,tissind)=mask(:)~=0;
            tiss_weight(:,tissind)=mask(:)/sum(mask(:));
        end;
        
        % Load up EPIs        
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
        
        % Don't collapse by tissue before saving timecourses
        % But do sort by tissue
        tiss_timecourse=cell(length(tissues),1);
        for tissind=1:length(tissues)
            tiss_timecourse{tissind}=Yepi_timecourse(tiss_mask(:,tissind)~=0,:);
        end;
        sesspth=aas_getsesspath(aap,subj,sess);
        tiss_timecourse_fn=fullfile(sesspth,'tiss_timecourse.mat');
        save(tiss_timecourse_fn,'tiss_timecourse','tiss_weight','tiss_mask');
        aap=aas_desc_outputs(aap,'session',[subj sess],'tiss_timecourse',tiss_timecourse_fn);
        
        %% Calculate power spectrum for each tissue compartment, weighted
        powerspect=abs(fft(Yepi_timecourse,[],2));
        powerspect=powerspect(:,1:size(powerspect,2)/2);
        tiss_powerspect=tiss_weight'*powerspect;
        
        nyquist=0.5/TR;
        freq=linspace(0,nyquist,size(tiss_powerspect,2));
        tiss_powerspect_fn=fullfile(sesspth,'tiss_powerspect.mat');
        save(tiss_powerspect_fn,'tiss_powerspect','freq');
        aap=aas_desc_outputs(aap,'session',[subj sess],'tiss_powerspect',tiss_powerspect_fn);
        
        %% Plot
        figure(1);
        plot(freq(2:end),tiss_powerspect(:,2:end));
        xlabel('Freq (Hz)');
        ylabel('Power');
        legend(tissues);
  
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;



