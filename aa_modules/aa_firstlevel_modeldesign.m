function [aap,resp]=aamod_firstlevel_modeldesign(aap,task,subj)

resp='';

switch task
    case 'report'

    case 'doit'
        % Movement regressors
        [moves, mnames] = aas_movPars(aap, subj, aap.tasklist.currenttask.settings.moveMat);


        % Other covariates ("compartment" regressors, globals, etc.)
        % [coming soon]



        % Get basis functions from task settings
        SPM.xBF = aap.tasklist.currenttask.settings.xBF;

        firstsess = aap.acq_details.selected_sessions(1);

        % if TR is manually specified, use that, otherwise try to get from DICOMs.
        if isfield(aap.tasklist.currenttask.settings,'TR') && ~isempty(aap.tasklist.currenttask.settings.TR)
            SPM.xY.RT = aap.tasklist.currenttask.settings.TR;
        else
            % Get TR from DICOM header checking they're the same for all sessions
            for sess = aap.acq_details.selected_sessions
                DICOMHEADERS=load(aas_getfiles_bystream(aap, subj, sess, 'epi_header'));
                try
                    TR=DICOMHEADERS.DICOMHEADERS{1}.volumeTR;
                catch
                    % [AVG] This is for backwards compatibility!
                    TR=DICOMHEADERS.DICOMHEADERS{1}.RepetitionTime/1000;
                end
                if (sess==firstsess)
                    SPM.xY.RT = TR;
                else
                    if (SPM.xY.RT~=TR)
                        aas_log(aap, true, sprintf('Session %d has different TR from earlier sessions, they can''t be in the same model.', sess));
                    end
                end
            end
        end

        % NB Previous versions tried to set T0 using sliceorder; I've left this out
        % assuming that users should either use the SPM defaults (which are in the XML
        % file) or set this manually.



        SPM.xGX.iGXcalc = aap.tasklist.currenttask.settings.global; % 'None';
        SPM.xVi.form = aap.tasklist.currenttask.settings.autocorrelation; % 'AR(1)';



    case 'checkrequirements'

    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end