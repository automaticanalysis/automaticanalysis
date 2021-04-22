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
        
        % Load up list of images and headers for each subject in this session
        V={};
        for subj=1:length(aap.acq_details.subjects)
            fn{subj}=aas_getimages_bystream(aap,subj,sess,'epi');
            % Load up headers and truncate if necessary
            [V{subj},~,~,thissubj_nimgs]=aas_spm_vol(fn{subj});             

            
            if aap.tasklist.currenttask.settings.truncatevolumes
                nimgs=min(nimgs,thissubj_nimgs);
            else
                if (subj>1 && nimgs~=thissubj_nimgs)
                    aas_log(aap,true,sprintf('Corresponding sessions of each subject must have the same number of EPI sessions - not true for session %s of subject %s as got %d not %d scans.',aap.acq_details.sessions(subj).name,aap.acq_details.subjects(subj).subjname,size(fn{subj},1),nimgs));
                end;                
            end;
        end;
        
        % Truncate headers if necessary
        for subj=1:length(aap.acq_details.subjects)
            V{subj}=V{subj}(1:nimgs);
        end;
        
        % Take average of each timepoint across subjects
        Ytot=[];
        Vout=V{1}(1);
        Vout.dt(1)=spm_type('float32');
        Vout.fname=fullfile(outdir,'meanepiacrosssubjects_4D.nii');  % Write file out a single 4D nii

        for k=1:nimgs
            for subj=1:length(aap.acq_details.subjects)
                if (subj==1)
                    Ytot=spm_read_vols(V{subj}(k));
                else
                    Ytot=Ytot+spm_read_vols(V{subj}(k));
                end;
            end;
            Ytot=Ytot/length(aap.acq_details.subjects);
            % k th image of 4D series
            Vout.n=[k 1];
            spm_write_vol(Vout,Ytot);
        end;
        

        aap=aas_desc_outputs(aap,'isc_session',sess,'meanepitimecourse',Vout.fname);
        
    otherwise
end;
