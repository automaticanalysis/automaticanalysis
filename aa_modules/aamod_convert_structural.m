function [aap,resp] = aamod_convert_structural(aap, task, subjInd)
% Generic m file for converting structural image (T1, T2).
%
% Used with aamod_convert_t1.xml, for example.
%
% Based on original aamod_copystructural.

resp='';

switch task       
    case 'description'
        resp='Structural DICOM to nifti';
        
    case 'summary'                
        
    case 'report'
        
    case 'doit'
    
        % Determine the type of image we are expecting (T1, T2, etc.) from
        % the input and output streams specified in the XML file.
        
        inStream = aap.tasklist.currenttask.inputstreams.stream{1};
        outStreamImg = aap.tasklist.currenttask.outputstreams.stream{1};
        outStreamDicom = aap.tasklist.currenttask.outputstreams.stream{2};
              
        
        % Convert the image and describe the output
        [aap convertedfns dcmhdr] = aas_convertseries_fromstream(aap, subjInd, inStream);                
        aap = aas_desc_outputs(aap, subjInd, outStreamImg, convertedfns);
        
        % save DICOM headers and describe the output
        subjpath = aas_getsubjpath(aap, subjInd);
        dcmhdrfn = fullfile(subjpath, sprintf('%s.mat', outStreamDicom));
        save(dcmhdrfn, 'dcmhdr');
        aap = aas_desc_outputs(aap, subjInd, outStreamDicom, dcmhdrfn);
 
    case 'checkrequirements'
        
    otherwise
        aas_log(aap, 1, sprintf('Unknown task %s', task));

end

