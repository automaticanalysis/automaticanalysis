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
        
        % It's possible that the sructural series contains more than one
        % volume, e.g., multiphase acquisition.  Maybe we only want to
        % carry a subset, or a single volume, forward.
        
        structVols = 1 : size(convertedfns, 2); % default to all volumes
        if isfield(aap.tasklist.currenttask.settings, 'struct_vols') && ~isempty(aap.tasklist.currenttask.settings.struct_vols)

            % Could be specified as a string, e.g., [1:4]
            if ischar(aap.tasklist.currenttask.settings.struct_vols)
                structVols = str2num(aap.tasklist.currenttask.settings.struct_vols);
            else
                structVols = aap.tasklist.currenttask.settings.struct_vols;
            end
            
            % If any of the specififed structural volumes don't exist in
            % the converted files, ignore this field.
            if any(~ismember(structVols, 1:size(convertedfns, 2)))
                structVols =  1:size(convertedfns, 2);
            end  
        end
        
        aap = aas_desc_outputs(aap, subjInd, outStreamImg, convertedfns(structVols, :));
        
        % save DICOM headers and describe the output
        subjpath = aas_getsubjpath(aap, subjInd);
        dcmhdrfn = fullfile(subjpath, sprintf('%s.mat', outStreamDicom));
        save(dcmhdrfn, 'dcmhdr');
        aap = aas_desc_outputs(aap, subjInd, outStreamDicom, dcmhdrfn);
 
    case 'checkrequirements'
        
    otherwise
        aas_log(aap, 1, sprintf('Unknown task %s', task));

end

