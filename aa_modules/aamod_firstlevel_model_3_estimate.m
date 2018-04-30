% AA module - first level statistics 3: estimate model (spm_spm)
%
% Split version of aamod_firstlevel_model - useful for inserting custom model
% options in between stages. Replaces aamod_firstlevel_model as a whole. Stages:
%
% aamod_firstlevel_model_1_config - add events, set up basic SPM struct
% aamod_firstlevel_model_2_convolve - generate convolved design matrix
% (spm_fmri_spm_ui)
% aamod_firstlevel_model_3_estimate - fit the model (spm_spm)

function [aap,resp]=aamod_firstlevel_model_3_estimate(aap,task,subj)

resp='';

switch task
    case 'report'

    case 'doit'
        % Get subject directory
        cwd=pwd;
        spmpath = aas_getfiles_bystream(aap,subj,'firstlevel_spm');
        load(spmpath);
        anadir = fileparts(spmpath);
        cd(anadir)
        % always a good idea to keep the swd current to avoid accidentally
        % overwriting across modules
        SPM.swd = anadir;
        
        %%%%%%%%%%%%%%%%%%%
        %% ESTIMATE MODEL%%
        %%%%%%%%%%%%%%%%%%%
        % avoid overwrite dialog
        prevmask = spm_select('List',SPM.swd,'^mask\..{3}$');
        if ~isempty(prevmask)
            for ind=1:size(prevmask,1)
                spm_unlink(fullfile(SPM.swd, prevmask(ind,:)));
            end;
        end
                
        SPMest = spm_spm(SPM);
        
        % Saving Residuals
        if isfield(aap.tasklist.currenttask.settings,'writeresiduals') && ~isempty(aap.tasklist.currenttask.settings.writeresiduals)
            aas_log(aap,false,'Writing residuals...');
            VRes = spm_write_residuals(SPMest,aap.tasklist.currenttask.settings.writeresiduals);
            for s = 1:numel(SPM.nscan)
                sesspath = aas_getsesspath(aap,subj,subjSessionI(s));
                fres = char({VRes(sum(SPM.nscan(1:s-1))+1:sum(SPM.nscan(1:s))).fname}');
                for f = 1:size(fres,1)
                    movefile(fres(f,:),sesspath);
                end
                residuals{subjSessionI(s)} = horzcat(repmat([sesspath filesep],[SPM.nscan(s) 1]),fres);
                if aap.options.NIFTI4D
                    spm_file_merge(residuals{subjSessionI(s)},fullfile(sesspath,sprintf('Res-%04d.nii',subjSessionI(s))),0,SPM.xY.RT);
                    fres = cellstr(residuals{subjSessionI(s)});
                    residuals{subjSessionI(s)} = fullfile(sesspath,sprintf('Res-%04d.nii',subjSessionI(s)));                
                    delete(fres{:});
                end                
            end
            aas_log(aap,false,'\bDone.');
        end        
        
        %% Describe outputs
        cd (cwd);

        % Describe outputs
        %  firstlevel_spm
        aap=aas_desc_outputs(aap,subj,'firstlevel_spm',fullfile(anadir,'SPM.mat'));

        %  firstlevel_betas (includes related statistical files)
        allbetas=dir(fullfile(anadir,'beta_*'));
        allbetas=vertcat(allbetas,...
            dir(fullfile(anadir,'ResMS.*')),...
            dir(fullfile(anadir,'RPV.*')));
        firstlevelmasking = ~all(isinf(SPM.xX.X));
        if firstlevelmasking
            allbetas=vertcat(allbetas,...
                dir(fullfile(anadir,'mask.*')));
            mask = dir(fullfile(anadir,'mask.*'));
            mask=strcat(repmat([anadir filesep],[numel(mask) 1]),char({mask.name}));
            aap=aas_desc_outputs(aap,subj,'firstlevel_brainmask',mask);            
        end
        betafns=strcat(repmat([anadir filesep],[numel(allbetas) 1]),char({allbetas.name}));
        aap=aas_desc_outputs(aap,subj,'firstlevel_betas',betafns);
        
        if isfield(aap.tasklist.currenttask.settings,'writeresiduals') && ~isempty(aap.tasklist.currenttask.settings.writeresiduals)
            for s = subjSessionI
                aap=aas_desc_outputs(aap,subj,s,'epi',residuals{s});
            end
        end

        %% DIAGNOSTICS...
        if ~isempty(SPM.xX.iC) % model is not empty
            h = firstlevelmodelStats(anadir, [], spm_select('FPList',anadir,'^mask.*'));
            print(h.regs,'-djpeg','-r150', fullfile(aas_getsubjpath(aap,subj), 'diagnostic_aamod_firstlevel_model_regs.jpg'));
            print(h.betas,'-djpeg','-r150', fullfile(aas_getsubjpath(aap,subj), 'diagnostic_aamod_firstlevel_model_betas.jpg'));
            
            close(h.regs)
            close(h.betas)
        end
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
