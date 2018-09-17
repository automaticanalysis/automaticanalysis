% AA module - collapse predictors in SPM according to the struct returned by
% calling collapsefun(SPM). 
%
% collapsestruct should be a struct array with one entry for each new regressor,
% and the fields name and targetnames (names of the regressors that should be
% collapsed over). Generally you would use this module would go between
% aamod_firstlevel_model_1_config and aamod_firstlevel_model_2_convolve in the
% task list.
%
% [aap,resp]=aamod_firstlevel_model_collapsepredictors(aap,task,subj)
function [aap,resp]=aamod_firstlevel_model_collapsepredictors(aap,task,subj)

resp='';

switch task
    case 'report'
    case 'doit'
        %get subject directory
        spmpath = aas_getfiles_bystream(aap,subj,'firstlevel_spm');
        ts = aap.tasklist.currenttask.settings;
        oldSPM = load(spmpath);
        oldSPM = oldSPM.SPM;
        % output struct
        SPM = oldSPM;
        SPM.Sess(:) = [];
        collapsestruct = feval(ts.collapsefun,oldSPM);
        for sess = 1:numel(oldSPM.Sess)
            % (don't) handle covariates
            if ~isempty(oldSPM.Sess(sess).C.C)
                aas_log(aap,1,'covariates are not currently supported');
            end
            SPM.Sess(sess).C = oldSPM.Sess(sess).C;
            % task regressors
            names = [oldSPM.Sess(sess).U.name];
            for thisu = collapsestruct(:)'
                newname = thisu.name;
                [~,targetind] = intersect(names,thisu.targetnames);
                if numel(targetind) ~= numel(thisu.targetnames)
                    aas_log(aap,1,'did not find the correct number of regressors');
                end
                newons = [oldSPM.Sess(sess).U(targetind).ons];
                newdur = [oldSPM.Sess(sess).U(targetind).dur];
                % (don't) handle modulators
                Ptest = [oldSPM.Sess(sess).U(targetind).P];
                if ~all(strcmp('none',{Ptest.name}))
                    aas_log(aap,1,'parametric modulators are not currently supported');
                end
                SPM.Sess(sess).U(end+1) = struct('ons',newons,...
                    'dur',newdur,'name',{{newname}},'P',Ptest(1),...
                    'orth',1);
            end
        end
        % these fields we will quietly ignore (usually set by
        % aamod_firstlevel_model_1_config)
        if isfield(SPM.xX,'iG')
            aas_log(aap,0,'SPM.xX.iG splitting is not supported - ignoring this field');
        end
        if isfield(SPM.xX,'iC')
            aas_log(aap,0,'SPM.xX.iC splitting is not supported - ignoring this field');
        end
        % Describe outputs
        %  firstlevel_spm
        save(spmpath,'SPM');
        aap=aas_desc_outputs(aap,subj,'firstlevel_spm',spmpath);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
