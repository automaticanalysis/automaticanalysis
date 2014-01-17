%aamod_fconn_computematrix Compute functional connectivity matrix
%
%  aamod_fc_computematrix(aap, subjInd, task) computes a functional
%  connectivity matrix based on timecourse relationships across seed regions.
%
%  settings.correlationtype can be set to 'pearson' (default), 'spearman',
%  or 'kendell'.

function [aap, resp] = aamod_fconnmatrix_seedseed(aap, task, subjInd)

resp='';

switch task
    case 'report'
        
    case 'doit'                        
        settings = aap.tasklist.currenttask.settings;                
        
        % Get image names
        %spmName = aas_getfiles_bystream(aap, subjInd, 'firstlevel_spm');

        
        % Loop through sessions and get vois for each
        
        sessions = aap.acq_details.selected_sessions;        
        nSessions = length(sessions);
                        
        for sessInd = 1:nSessions
            
            thisSess = sessions(sessInd);            
            sessDir = aas_getsesspath(aap, subjInd, thisSess);
            
            voiFile = aas_getfiles_bystream(aap, subjInd, thisSess, 'voi'); % timecourses
            load(voiFile)
            
            % Make sure all sessions have the same number of regions
            if sessInd==1
                nVoi = length(vois);
            else
                if length(vois)~=nVoi
                    aas_log(aap, true, sprintf('Different number of VOIs in session %i. Expecting same VOIs in each session.', thisSess));
                end
            end
            
            nTime = length(vois(1).data.mean); % how many time points
            dataMean{sessInd} = zeros(nTime, nVoi);
            dataEigen{sessInd} = zeros(nTime, nVoi);
            
            for voiInd = 1:nVoi
                dataMean{sessInd}(:,voiInd) = vois(voiInd).data.mean;
                dataEigen{sessInd}(:,voiInd) = vois(voiInd).data.firsteigenvalue;                
            end                        
        end
        
        
        %TODO: ignore spikes before performing correlation -JP
                                        
        % Do correlations across all volumes (concatenate) or split by
        % session and then average?
        if settings.concatenate
            [rMean, pMean] = corr(vertcat(dataMean{:}), 'type', settings.correlationtype);
            [rEigen, pEigen] = corr(vertcat(dataEigen{:}), 'type', settings.correlationtype);
        else
           aas_log(aap, true, 'Not concatenating is not implemented for now, sorry'); 
            
        end                        
        
        % output
        %[pth, nm, ext] = fileparts(spmName);
        pth = aas_getsubjpath(aap, subjInd);
        subName = aap.acq_details.subjects(subjInd).mriname;
        fconnName = sprintf('%s-fconn%s.mat', subName, settings.matrixsuffix);
        fconnPath = fullfile(pth, fconnName);
        
        % output matrix
        fconn = [];
        fconn.subjectName = subName;
        fconn.regionNames = {vois.name};
        fconn.rMean = rMean;
        fconn.pMean = pMean;
        fconn.rMeanZ = .5 * log((1+rMean)./(1-rMean)); % Fisher transform: https://en.wikipedia.org/wiki/Fisher_transformation
        fconn.rEigen = rEigen;
        fconn.pEigen = pEigen;
        fconn.rEigenZ = .5 * log((1+rEigen)./(1-rEigen)); % Fisher transform
        
        fconn.concatenatedbysession = settings.concatenate;
        fconn.correlationtype = settings.correlationtype;                
        fconn.README = 'First iteration of connectivity matrix. What''s the greek letter before alpha?';
        
        save(fconnPath, 'fconn');
        
        % If there's a matrix farm, make a link
        if ~isempty(settings.matrixfarm)
            try
                if ~isdir(settings.matrixfarm)
                    mkdir(settings.matrixfarm);
                end
                
                farmPath = fullfile(settings.matrixfarm, fconnName);
                
                % Copy or link, depending
                if settings.matrixfarmlink==0
                    system(sprintf('cp %s %s', fconnPath, farmPath));
                else
                    system(sprintf('ln -s %s %s', fconnPath, farmPath));
                end
            catch
                aas_log(aap, false, sprintf('Tried to link .mat file to matrixfarm %s but it did not work.', settings.matrixfarm));
            end
        end
                
        % Describe output
        aap = aas_desc_outputs(aap, subjInd, 'firstlevel_fconn_matrix_mat', fconnPath);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end




