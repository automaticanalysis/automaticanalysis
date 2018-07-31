% AA module - across subject movie correlations
% Rhodri Cusack MRC CBU Cambridge Dec 2010
%
% At each voxel, correlate an individual subject's timecourse with the mean
% timecourse across all of the other subjects
%

function [aap,resp]=aamod_moviecorr_meantimecourse(aap,task,sess,subj)

resp='';

switch task
    case 'report'
        
    case 'doit'
        nsub=length(aap.acq_details.subjects);
        subj_imgs=aas_getimages_bystream(aap,subj,sess,'epi');
        % This is a subject-level thing
        mean_imgs=aas_getfiles_bystream(aap,'isc_session',sess,'meanepitimecourse');
        % Truncate
        if (~aap.tasklist.currenttask.settings.truncatevolumes && size(subj_imgs,1)~=size(mean_imgs,1))
            aas_log(aap,true,sprintf('Expected same number of images in epi and meanepitimecourse\n'));
            Vfirstmean=spm_vol(mean_imgs(1,:));
            Vfirstsubj=spm_vol(subj_imgs(1,:));
            if (any(Vfirstmean.dim~=Vfirstsubj.dim) || any(Vfirstmean.mat~=Vfirstsubj.mat))
                aas_log(aap,true,sprintf('Need same image dimensions and space for mean and individual subject.'));
            end;
        end;
        
        % This will return a structure  from 4D or 3D inputs
        V1 =aas_spm_vol(subj_imgs);
        V2 =aas_spm_vol(mean_imgs);

        nimg=length(V1);
        nsubj=length(aap.acq_details.subjects);

        for imgind=1:nimg
            Y1=spm_read_vols(V1(imgind));
            Y2=spm_read_vols(V2(imgind));
            
            % Subtract this subject's contribution, if the mean includes
            % them
            if (aap.tasklist.currenttask.settings.meanincludeseachsubject) 
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
        Vout=V1(1);
        Vout.fname=fullfile(aas_getpath_bydomain(aap,'isc_subject',[sess subj]),'moviecorr_loo.nii');
        Vout.dt(1)=spm_type('float32');
        spm_write_vol(Vout,corr(:,:,:));
        aap = aas_desc_outputs(aap,'isc_subject',[sess subj],'moviecorr_loo',Vout.fname);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;














