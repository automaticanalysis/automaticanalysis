% automatic analysis module -
% Calculates the mean of the EPIs for a single session
% This module has domain isc_session, so unusually, sess is the only
% parameter
% Rhodri Cusack 2011-2015, BMI Western, Canada

function [aap resp]=aamod_meanepitimecourse(aap,task,sess)
resp='';
switch task
    case 'doit'
        
        outputs={};
        nimgs=inf;
        % Output directory is an session directory at the study
        % level (unusual)
        outdir=aas_getpath_bydomain(aap,'isc_session',sess);
        aap=aas_makedir(aap,outdir);
        
        % Load up list of images for each subject in this session
        for i=1:length(aap.acq_details.subjects)
            fn{i}=aas_getimages_bystream(aap,i,sess,'epi');
            if aap.tasklist.currenttask.settings.truncatevolumes
                nimgs=min(nimgs,size(fn{i},1));
            else
                if (i>1 && nimgs~=size(fn{i},1))
                    aas_log(aap,true,sprintf('Corresponding sessions of each subject must have the same number of EPI sessions - not true for session %s of subject %s as got %d not %d scans.',aap.acq_details.sessions(i).name,aap.acq_details.subjects(i).mriname,size(fn{i},1),nimgs));
                end;
                nimgs=size(fn{i},1);
            end;
        end;
        % Take average of each timepoint across subjects
        for k=1:nimgs
            for i=1:length(aap.acq_details.subjects)
                V=spm_vol(fn{i}(k,:));
                if (i==1)
                    Ytot=spm_read_vols(V);
                else
                    Ytot=Ytot+spm_read_vols(V);
                end;
            end;
            
            Ytot=Ytot/length(aap.acq_details.subjects);
            V.fname=fullfile(outdir,sprintf('meanepiacrosssubjects_%05d.nii',k));
            V.dt(1)=spm_type('float32');
            spm_write_vol(V,Ytot);
            outputs=[outputs V.fname];
        end;
        aap=aas_desc_outputs(aap,'isc_session',sess,'meanepitimecourse',outputs);
        
    otherwise
end;
