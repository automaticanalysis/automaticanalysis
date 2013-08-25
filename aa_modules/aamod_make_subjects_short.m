% AA initialisation module - make subject filenames for structural image
% format controlled by aap.directory_conventsion.subject_filenames_format
%        1: Make short subject filenames (e.g., CBU030511)
%        2: Make ordinal short subject filenames (e.g., S01, S02...)'
%        0: User specified short subject filenames';

% Rhodri Cusack MRC CBU Cambridge 2004


function [aap,resp]=aamod_make_subjects_short(aap,task)

resp='';

switch task
    case 'domain'
        resp='study';   % this module needs to be run once per study
    case {'description','summary'}
        switch(aap.directory_conventions.subject_filenames_format)
            case 1
                resp='Make short subject filenames (e.g., CBU030511)';
            case 2
                resp='Make ordinal short subject filenames (e.g., S01, S02...)';
            case 0
                resp='User specified short subject filenames';
            case 3
                resp='Make short subject filenames (e.g., CBU030511 unless already specified by user';
             case 4
                resp='Just filter for any nifti files';

            otherwise
                aas_log(1,sprintf('Unknown subject directory format (aap.directory_conventions.subjectfn_format=%d)',aap.directory_conventions.subjectfn_format));
        end;
    case 'doit'
        switch(aap.directory_conventions.subject_filenames_format)
            case {1,3}

                for i=1:length(aap.acq_details.subjects)
                    if (aap.directory_conventions.subject_filenames_format==1 || (aap.directory_conventions.subject_filenames_format==3 && isempty(aap.acq_details.subjects(i).structuralfn)))

                        % djm: added 1st two options to cope with fully
                        % specified or missing mrinames (160708)
                        if strcmpi('missing',aap.acq_details.subjects(i).mriname)
                            continue
                        elseif exist(aap.acq_details.subjects(i).mriname,'file')
                            [junk aap.acq_details.subjects(i).structuralfn]=fileparts(aap.acq_details.subjects(i).mriname);
                        else
                            % other sites will probably need to change
                            % these
                            % parses things like TESTCBU060500ANYTHING into CBU060500
                            tmp=strtok(aap.acq_details.subjects(i).mriname,'0123456789');  % find part before first digit (that'll be 'TESTCBU' in above example)
                            pos=length(tmp)-2;                                          % two characters less than length (giving 7-2=5 in example)
                            tmp=aap.acq_details.subjects(i).mriname(pos:(pos+8)); % pick out 9 characters from this (characters 5 to 13 'CBU060500' in example)

                            aap.acq_details.subjects(i).structuralfn=tmp;
                        end
                    end;
                end;
            case 2
                for i=1:length(aap.acq_details.subjects)
                    aap.acq_details.subjects(i).structuralfn=sprintf('S%2d',i);
                end;
            case 4
                for i=1:length(aap.acq_details.subjects)
                    aap.acq_details.subjects(i).structuralfn='*';
                end;
        end;


    case 'checkrequirements'
        for i=1:length(aap.acq_details.subjects)
            if isempty(aap.acq_details.subjects(i).mriname)
                aas_log(aap,1,sprintf('No filename specified for subject %d\n',i));
            end;
        end;

    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
end;
