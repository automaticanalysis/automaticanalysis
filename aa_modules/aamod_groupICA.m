% Run a group ICA analysis with 'gift'

function [aap, resp] = aamod_groupICA(aap, task)

resp='';

switch task
    case 'domain'
        resp='study';   % this module needs to be run once per study
        
    case 'description'
        resp='Group ICA analysis';
        
    case 'summary'
        subjpath = aas_getsubjpath(subj);
        resp = sprintf('Gift model %s\n', subjpath);
        
    case 'report'
        
    case 'doit'

        files = struct('name', {});
        
        for sub = 1 : length(aap.acq_details.subjects)
            
            [numSess, sessInds] = aas_getN_bydomain(aap, 'session', sub);
            
            for ses = sessInds
                files(end+1) = struct('name', aas_getimages_bystream(aap, sub, ses, 'epi'));
            end

        end
         
        % Prepare variables expected by GIFT
        modalityType = 'fMRI';
        numOfSess = length(aap.acq_details.selected_sessions);
        numOfSub = length(aap.acq_details.subjects);
        SPMFiles = struct('name', []);
        
        studyPath = aas_getstudypath(aap);
        dataMatFile = fullfile(studyPath, 'Subject.mat');
        save(dataMatFile, 'files', 'modalityType', 'numOfSess', 'numOfSub', 'SPMFiles');
        aap = aas_desc_outputs(aap, 'gift', dataMatFile);
        

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end

end