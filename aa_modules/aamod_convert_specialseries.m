% AA module - Converts special series to NIFTI format
% Rhodri Cusack MRC CBU Cambridge Nov 2005

function [aap,resp]=aamod_convert_specialseries(aap,task,subj,sess)

resp='';

switch task
    case 'report'
    case 'doit'
        for inpstream = aas_getstreams(aap,'input')            
            streamname = strrep(inpstream{1},'dicom_','');
            
            [aap, convertedfns, dcmhdr] = aas_convertseries_fromstream(aap, subj, sess, ['dicom_' streamname]);
            
            outstream = {};
            % Restructure outputs!
            if iscell(convertedfns{1}) % subdirs
                for c = 1:length(convertedfns)
                    outstream = [outstream; convertedfns{c}];
                    dcmhdr{c} = dcmhdr{c}{1}; 
                end
            else
                outstream = convertedfns;
            end
            
            aap = aas_desc_outputs(aap, 'special_session', [subj, sess], streamname, outstream);
            dcmhdrfn = fullfile(aas_getsesspath(aap,subj,sess),[streamname '_dicom_headers.mat']);
            save(dcmhdrfn,'dcmhdr');
            aap = aas_desc_outputs(aap, 'special_session', [subj, sess], [streamname '_dicom_header'], dcmhdrfn);
        end
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
