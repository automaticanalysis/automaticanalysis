% AA module - split first-level SPM module sessions into sub-runs. Should
% typically be placed in between aamod_firstlevel_model_1_config and
% aamod_firstlevel_model_2_convolve in the task list.

function [aap,resp]=aamod_firstlevel_model_subrunsplit(aap,task,subj)

resp='';

switch task
    case 'report'
    case 'doit'

        spmpath = aas_getfiles_bystream(aap,subj,'firstlevel_spm');
        load(spmpath);
        anadir = fileparts(spmpath);
        % always a good idea to keep the swd current to avoid accidentally
        % overwriting across modules
        SPM.swd = anadir;

        split = aap.tasklist.currenttask.settings.split;

        % onsets. We need to subtract an offset value to make the new
        % first volume be time 0
        % (- -1 because we go from 1-based indexing to time)
        % scan timing
        ind2onset = @(x)x-1;
        if strcmp(SPM.xBF.UNITS,'secs')
            % second timing
            ind2onset = @(x)(x-1) * SPM.xY.RT;
        end

        newSPM = SPM;
        newSPM.xY.P = [];
        newSPM.nscan = [];
        newSPM.Sess = struct;
        scanoff = 0;
        newsess = 0;
        for sess = 1:numel(SPM.Sess)
            % volumes in this session
            scanind = scanoff + (1:SPM.nscan(sess));
            % update running index
            scanoff = scanind(end);
            for thissplitc = split(:)'
                thissplit = thissplitc{1};
                newsess = newsess + 1;
                % data
                splitind = scanind(thissplit(1):thissplit(2));
                newSPM.nscan(newsess) = numel(splitind);
                newSPM.xY.P = [newSPM.xY.P; SPM.xY.P(splitind,:)];
                % covariates
                if isfield(SPM.Sess(sess),'C')
                    newSPM.Sess(newsess).C = SPM.Sess(sess).C;
                    if ~isempty(SPM.Sess(sess).C.C)
                        newSPM.Sess(newsess).C.C = newSPM.Sess(newsess).C(thissplit(1):thissplit(2),:);
                    end
                end
                % shift all onsets to account for new time zero 
                onsoffset = ind2onset(thissplit(1));
                % and trim off onsets that are now outside the run
                maxtime = ind2onset(thissplit(2));
                for conind = 1:numel(SPM.Sess(sess).U);
                    % shift start time
                    thisu = SPM.Sess(sess).U(conind);
                    thisu.ons = thisu.ons - onsoffset;
                    badtrials = thisu.ons < 0 | thisu.ons > maxtime;
                    thisu.ons(badtrials) = [];
                    thisu.dur(badtrials) = [];
                    if ~all(strcmp(thisu.P.name,'none'))
                        aas_log(aap,true,['parametric modulators (P field ' ...
                        'in SPM.Sess.U) are not currently supported']);
                    end
                    newSPM.Sess(newsess).U(conind)= thisu;
                end
                if isfield(SPM,'xX') && isfield(SPM.xX,'K')
                    newSPM.xX.K(newsess) = SPM.xX.K(sess);
                end
            end
        end
        if isfield(SPM,'xX')
            % these fields we will quietly ignore (usually set by
            % aamod_firstlevel_model_1_config)
            if isfield(SPM.xX,'iG')
                aas_log(aap,0,'SPM.xX.iG splitting is not supported - ignoring this field');
            end
            if isfield(SPM.xX,'iC')
                aas_log(aap,0,'SPM.xX.iC splitting is not supported - ignoring this field');
            end
        end

        %oldSPM = SPM;
        SPM = newSPM;
        save(spmpath,'SPM');
        aap=aas_desc_outputs(aap,subj,'firstlevel_spm',spmpath);
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end
