% AA module - Converts special series to NIFTI format
% Rhodri Cusack MRC CBU Cambridge Nov 2005

function [aap,resp]=aamod_convert_specialseries(aap,task,subj)

resp='';

switch task
    case 'report'
    case 'doit'
        subjpath=aas_getsubjpath(aap,subj);

        streamname = strrep(aas_getstreams(aap,'output'),'dicom_',''); streamname = streamname{1};

        [aap, convertedfns, dcmhdr] = aas_convertseries_fromstream(aap, subj, ['dicom_' streamname]); 
        
        outstream = {};
        % Restructure outputs!
        for c = 1:length(convertedfns)
            outstream = [outstream; convertedfns{c}];
        end
        
        % reorder if there is a baseline session and it is not the first (for coregistration)
        ind = cell_index(outstream,'baseline'); ind = ind(1);
        if ind && (ind ~=1)
            outstream = outstream([ind 1:ind-1 ind+1:end]);
        end
        
        aap = aas_desc_outputs(aap, 'subject', subj, streamname, outstream);
        dcmhdrfn = fullfile(subjpath,'dicom_headers.mat');
        save(dcmhdrfn,'dcmhdr');
        aap = aas_desc_outputs(aap, 'subject', subj, [streamname '_dicom_header'], dcmhdrfn);
        
    case 'checkrequirements'

    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
