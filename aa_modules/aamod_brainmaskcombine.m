% AA module - Create a binary brainmask for a group of subjects based on
% individually-defined binary masks
%
% summer 2023 [MSJ] - modifed to use renamble instream
%

function [aap, resp] = aamod_brainmaskcombine(aap, task)
resp='';

switch task
    
    case 'report'
        
    case 'doit'
              
        mask_inputstream_struct = aap.tasklist.currenttask.inputstreams(1).stream{1};

        nSubj = length(aap.acq_details.subjects);
        
        for subjInd = 1:nSubj
            singleMasks{subjInd} = aas_getfiles_bystream(aap, subjInd, mask_inputstream_struct.CONTENT);         
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
