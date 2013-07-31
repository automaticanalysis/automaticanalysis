% AA module - Create a binary brainmask for a group of subjects based on
% individually-defined binary masks
%
function [aap, resp] = aamod_brainmaskcombine(aap, task)
resp='';

switch task
    case 'report'
        
    case 'doit'
        
        % how many subjects?
        nSubj = length(aap.acq_details.subjects);
        
        % for each subject, get their brain mask
        for subjInd = 1:nSubj
            singleMasks{subjInd} = aas_getfiles_bystream(aap, subjInd, 'brainmask');
        end
        
       Vi = spm_vol(strvcat(singleMasks));
       Vo = Vi(1);
       
       studyPath = aas_getstudypath(aap);
       outName = fullfile(studyPath, 'groupbinarybrainmask.nii');
       Vo.fname = outName;
       
       f = 'all(X)'; % all 1s
       flags = {1,0,0};
       
        
       Vo = spm_imcalc(Vi,Vo,f,flags);
       
       aap = aas_desc_outputs(aap, 'groupbrainmask', outName);
       
        
    case 'checkrequirements'
        
end
