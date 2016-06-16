% AA module - Converts special series to NIFTI format
% Rhodri Cusack MRC CBU Cambridge Nov 2005

function [aap,resp]=aamod_convert_specialseries(aap,task,subj,sess)

resp='';

switch task
    case 'report'
    case 'doit'

        streamname = strrep(aas_getstreams(aap,'output'),'dicom_',''); streamname = streamname{1};

        [aap, convertedfns, dcmhdr] = aas_convertseries_fromstream(aap, subj, sess, ['dicom_' streamname]); 
        
        outstream = {};
        % Restructure outputs!
        for c = 1:length(convertedfns)
            outstream = [outstream; convertedfns{c}];
        end
        
        aap = aas_desc_outputs(aap, 'special_session', [subj, sess], streamname, outstream);
        dcmhdrfn = fullfile(aas_getsesspath(aap,subj,sess),'dicom_headers.mat');
        save(dcmhdrfn,'dcmhdr');
        aap = aas_desc_outputs(aap, 'special_session', [subj, sess], [streamname '_dicom_header'], dcmhdrfn);
        
    case 'checkrequirements'

    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
