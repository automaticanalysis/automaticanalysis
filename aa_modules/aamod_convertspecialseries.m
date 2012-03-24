% AA module - Converts special series to NIFTI format
% Rhodri Cusack MRC CBU Cambridge Nov 2005

function [aap,resp]=aamod_convertspecialseries(aap,task,i)

resp='';

switch task
    case 'domain'
        resp='subject';   % this module needs to be run once per subject

    case 'description'
        resp='Convert special series';

    case 'summary'
        somethingtoconvert=false;
        if (length(aap.acq_details.specialseries)~=0)
            if (length(aap.acq_details.specialseries{i})==0)
                somethingtoconvert=true;
            end;
        end;
        if (somethingtoconvert)
            resp=sprintf('No special series for subject %s\n',aap.acq_details.subjects(i).mriname);
        else
            resp=sprintf('Converted special series for subject %s to %s\n', aap.acq_details.subjects(i).mriname,aap.directory_conventions.centralstore_structurals);
        end;

    case 'report'
    case 'doit'
        subjpath=aas_getsubjpath(aap,i);
        specialseriesdir=fullfile(subjpath,aap.directory_conventions.specialseriesdirname);
        aas_makedir(aap,specialseriesdir);

        for k=1:length(aap.acq_details.specialseries{i})
            subjpath=aas_getsubjpath(aap,i);
            dicomdirsearchpath=fullfile(aap.directory_conventions.rawdatadir,aap.acq_details.subjects(i).mriname,sprintf('Series_%03d*',aap.acq_details.specialseries{i}(k)));
            [aap dicomdatadir]=aas_findfiles(aap,dicomdirsearchpath,1,1);
            specialseriesdir_sub=fullfile(specialseriesdir,sprintf('series_%d',aap.acq_details.specialseries{i}(k)));
            aas_makedir(aap,specialseriesdir_sub);
            [aap]=aas_convertseries(aap,i,aap.acq_details.specialseries{i}(k),specialseriesdir_sub);
        end;
    case 'checkrequirements'

    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
