% AA module
% Performs .zil to .txt conversion of all DMDX data files
% You will also need the "getzil.awk" script that I wrote
% Matt Davis MRC CBU Cambridge June 2006

function [aap,resp]=aamod_convert_dmdx(aap,task,i,j)

resp='';

switch task
    case 'domain'
        resp='session';   % this module needs to be run once per session
        
    case 'description'
        resp=sprintf('Performing DMDX to text conversion of event files');
        
    case 'summary'
        resp=sprintf('Made directory %s and did dicom to nifti conversion of EPIs\n');
    case 'report'
    case 'doit'
        subjpath=aas_getsubjpath(aap,i);
        eventsdir=fullfile(subjpath,aap.directory_conventions.eventsdirname);
        
        % If doesn't exist then make events directory, as this is used to store done_flags        
        % but this module is called once for each session... problem?
        if (~length(dir(eventsdir)))
            [s w]=aas_shell(['mkdir ' eventsdir]);
            if (s)
                aas_log(aap,1,sprintf('Problem creating directory for session\n%s',sesspath));
            end;
        end;
        

        [aap]=aas_convertdmdx(aap,i,aap.acq_details.brukersessionnums{i}(j),eventsdir);
        
    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;


%_______________________________________________________________________

