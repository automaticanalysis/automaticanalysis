% AA module - Converts general series to NIFTI format
% Rhodri Cusack MRC CBU Cambridge Nov 2005

function [aap,resp]=aamod_convertphasecontrast(aap,task,i)

resp='';

switch task
    case 'domain'
        resp='subject';   % this module needs to be run once per subject

    case 'description'
        resp='Convert fieldmaps';

    case 'summary'
        if (length(aap.acq_details.siemensfieldmap{i})==0)
            resp=sprintf('No fieldmaps for subject %s\n',aap.acq_details.subjects(i).mriname);
        else
            resp=sprintf('Converted fieldmaps for subject %s to %s\n', aap.acq_details.subjects(i).mriname,aap.directory_conventions.centralstore_structurals);
        end;

    case 'report'
    case 'doit'
        subjpath=aas_getsubjpath(aap,i);
        fieldmapsdir=fullfile(subjpath,aap.directory_conventions.fieldmapsdirname);
        if (~length(dir(fieldmapsdir)))
            [s w]=aas_shell(['mkdir ' fieldmapsdir]);
            if (s)
                aas_log(aap,1,sprintf('Problem making directory%s',fieldmapsdir));
            end;
        end;
        fieldmapsmagdir=fullfile(fieldmapsdir,'rawmag');
        if (~length(dir(fieldmapsmagdir)))
            [s w]=aas_shell(['mkdir ' fieldmapsmagdir]);
            if (s)
                aas_log(aap,1,sprintf('Problem making directory%s',fieldmapsmagdir));
            end;
        end;
        fieldmapsphasedir=fullfile(fieldmapsdir,'rawphase');
        if (~length(dir(fieldmapsphasedir)))
            [s w]=aas_shell(['mkdir ' fieldmapsphasedir]);
            if (s)
                aas_log(aap,1,sprintf('Problem making directory%s',fieldmapsphasedir));
            end;
        end;

        if (length(aap.acq_details.siemensfieldmap{i})>2)
            aas_log(aap,1,'More than a pair of fieldmap series. You must choose one pair, and then specified the others as ignored as a final parameter of aas_addsubject command in your user script.\n');
        end;

        if (length(aap.acq_details.siemensfieldmap{i})==2)
            [aap]=aas_convertseries(aap,i,aap.acq_details.siemensfieldmap{i}(1),fieldmapsmagdir);
            [aap]=aas_convertseries(aap,i,aap.acq_details.siemensfieldmap{i}(2),fieldmapsphasedir);
        end;
    case 'checkrequirements'

    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
