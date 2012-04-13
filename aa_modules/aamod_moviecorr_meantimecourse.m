% AA module - across subject movie correlations
% Rhodri Cusack MRC CBU Cambridge Dec 2010
%
% At each voxel, correlate an individual subject's timecourse with the mean
% timecourse across all of the other subjects
%

function [aap,resp]=aamod_moviecorr_meantimecourse(aap,task,i,j)

resp='';

switch task
    case 'report'
        
    case 'doit'
        nsub=length(aap.acq_details.subjects);
        subj_imgs=aas_getimages_bystream(aap,i,j,'epi');
        % This is a subject-level thing
        mean_imgs=aas_getfiles_bystream(aap,'meanepitimecourse');
        % A bit messy - now filter for the session files in the study directory we care about
        mean_imgs=mean_imgs(strmatch(fullfile(aas_getstudypath(aap),aap.acq_details.sessions(j).name),mean_imgs),:);
        if (size(subj_imgs,1)~=size(mean_imgs,1))
            aas_log(aap,true,sprintf('Expected same number of images in epi and meanepitimecourse\n'));
            Vfirstmean=spm_vol(mean_imgs(1,:));
            Vfirstsubj=spm_vol(subj_imgs(1,:));
            if (any(Vfirstmean.dim~=Vfirstsubj.dim) || any(Vfirstmean.mat~=Vfirstsubj.mat))
                aas_log(aap,true,sprintf('Need same image dimensions and space for mean and individual subject.'));
            end;
        end;
        nimg=size(subj_imgs,1);
        nsubj=length(aap.acq_details.subjects);
        
        for imgind=1:nimg
            V1=spm_vol(subj_imgs(imgind,:));
            Y1=spm_read_vols(V1);
            V2=spm_vol(mean_imgs(imgind,:));
            Y2=spm_read_vols(V2);
            
            % Subtract this subject's contribution, if the mean includes
            % them
            if (aap.tasklist.currenttask.meanincludesthissubject)
                Y2=Y2-Y1/nsubj;
            end;
            
            % Stats
            if (imgind==1)
                Y1tot=Y1;
                Y1totsq=Y1.^2;
                Y2tot=Y2;
                Y2totsq=Y2.^2;
                Y1Y2tot=Y1.*Y2;
            else
                Y1tot=Y1tot+Y1;
                Y1totsq=Y1totsq+Y1.^2;
                Y2tot=Y2tot+Y2;
                Y2totsq=Y2totsq+Y2.^2;
                Y1Y2tot=Y1Y2tot+Y1.*Y2;
            end;
        end;
        corr=(nimg*Y1Y2tot-Y1tot.*Y2tot)./(sqrt(nimg*Y1totsq-Y1tot.^2).*sqrt(nimg*Y2totsq-Y2tot.^2));
        % Any correlations outside -1 to 1 must be division errors
        corr(corr<-1)=-1;
        corr(corr>1)=1;
        V1.fname=fullfile(aas_getsesspath(aap,i,j),'moviecorr_meantimecourse.nii');
        V1.dt(1)=spm_type('float32');
        spm_write_vol(V1,corr(:,:,:));
        aap = aas_desc_outputs(aap,i,j,'moviecorr_meantimecourse','moviecorr_meantimecourse.nii');
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;














