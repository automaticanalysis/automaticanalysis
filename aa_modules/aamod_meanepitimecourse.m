function [aap resp]=aamod_meanepitimecourse(aap,task)
resp='';
switch task
    case 'doit'
        % Unsually, sessions is the outer loop
        for j=aap.acq_details.selected_sessions
            outputs={};
            % Load up list of images for each subject in this session
            for i=1:length(aap.acq_details.subjects)
                fn{i}=aas_getimages_bystream(aap,i,j,'epi');
                if (i>1 && nimgs~=size(fn{i},1))
                    aas_log(aap,true,sprintf('Corresponding sessions of each subject must have the same number of EPI sessions - not true for session %s of subject %s as got %d not %d scans.',aap.acq_details.sessions(i).name,aap.acq_details.subjects(i).name,size(fn{i},1),nimgs));
                end;
                nimgs=size(fn{i},1);
            end;
            % Output directory is an session directory at the study
            % level (unusual)
            outdir=fullfile(aas_getstudypath(aap),aap.acq_details.sessions(j).name);
            aap=aas_makedir(aap,outdir);
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
            aap=aas_desc_outputs(aap,'meanepitimecourse',outputs);
        end;
    otherwise
end;
