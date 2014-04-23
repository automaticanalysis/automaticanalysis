% Apply bedpostx
% BEDPOSTX stands for Bayesian Estimation of Diffusion Parameters Obtained 
% using Sampling Techniques. The X stands for modelling Crossing Fibres. 
% bedpostx runs Markov Chain Monte Carlo sampling to build up distributions
% on diffusion parameters at each voxel. It creates all the files necessary
% for running probabilistic tractography

function [aap resp]=aamod_diffusion_bedpostx(aap,task)
global aaworker
resp='';


switch task
    case 'report'
    case 'doit'
        % Launch all of the bedpost jobs
        sessind=1;
        fslext=aas_getfslext(aap);
        bedpostpath={};
        diffusion_pth={};
        dsesspth={};
        
        for subjind=1:length(aap.acq_details.subjects)
            dsesspth{subjind}= aas_getpath_bydomain(aap,'diffusion_session',[subjind sessind]);
            bedpostpath{subjind}=[dsesspth{subjind} '.bedpostX'];
            
            % Change name to data.nii.gz for bedpost
            diffusion_pth{subjind}=fullfile(dsesspth{subjind}, ['data' fslext]);
            try
                movefile(aas_getfiles_bystream(aap,'diffusion_session',[subjind sessind],'diffusion_data'),diffusion_pth{subjind});
            catch
                aas_log(aap,false,sprintf('Couldn''t move, so assuming already moved.'));
            end;
            % Remove bet_ from start of bet_nodif_brain_mask etc as bedpost
            % doesn't like this
            betfiles=aas_getfiles_bystream(aap,'diffusion_session',[subjind sessind],'BETmask');
            for betind=1:size(betfiles,1)
                [pth nme ext]=aas_fileparts(betfiles(betind,:));
                if strcmp(nme(1:4),'bet_')
                    movefile(deblank(betfiles(betind,:)),deblank(fullfile(pth,[nme(5:end) ext])));
                end;
            end;
            
            % Don't redo if bedpost has completed already
            [aap,alldone]=aas_bedpostx_monitor(aap,diffusion_pth{subjind},bedpostpath{subjind},false,false);
            if alldone
                aas_log(aap,false,sprintf('Bedpost is complete, not rerunning: %s',bedpostpath{subjind}));
            else
                % Remove bedpost output path
                cmd=sprintf('rm -rf %s',bedpostpath{subjind});
                [s w]=aas_shell(cmd,true);
                
                cmd=sprintf('bedpostx %s --nf=2 --fudge=1  --bi=1000 --model=2',dsesspth{subjind})
                [s w]=aas_runfslcommand(aap,cmd);
                if s
                    aas_log(aap,true,sprintf('Error executing\n  %s\nof\n%s',cmd,w));
                end;
            end;
        end;
        
        % Now wait for all of these to finish
        isfirsttime=false; %true;
        subj_alldone=false(subjind);
        while any(~subj_alldone)
            for subjind=1:length(aap.acq_details.subjects)
                % And wait for it to finish!!
                if ~subj_alldone(subjind)
                    [aap,alldone]=aas_bedpostx_monitor(aap,diffusion_pth{subjind},bedpostpath{subjind},false,isfirsttime);
                    if alldone
                        % Now describe outputs
                        outstreams=aap.tasklist.currenttask.outputstreams;
                        for outind=1:length(outstreams.stream)
                            aap=aas_desc_outputs(aap,'diffusion_session_bedpostx',[subjind sessind],outstreams.stream{outind},fullfile(bedpostpath{subjind},[outstreams.stream{outind} fslext]));
                        end;
                        subj_alldone(subjind)=alldone;
                    end;
                end;
            end
        end
end
end

