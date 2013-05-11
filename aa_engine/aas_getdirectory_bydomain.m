% Automatic analysis - this file determines the names of each of the
% directory levels
%   domain='subject','session' etc
%   index= number of item

function [directory]=aas_getdirectory_bydomain(aap,domain,index)

switch (domain)
    case 'searchlight_package'
        directory=sprintf('searchlight_package_%d',index);
    
    case 'hyperalignment_searchlight_package'
        directory=sprintf(['hyperalignment_searchlight_packages' filesep '%d'],index);
        
    case {'splitsession_cv_fold','splitsession_cv_fold_hyper'}
        directory=sprintf('splitsession_cv_fold_%d',index);
        
    case 'session'
        directory=aap.acq_details.sessions(index).name;
        
        
    case {'subject','hyperalignment_subject'}
        switch (aap.directory_conventions.subject_directory_format)
            case 1
                if (isfield(aap.acq_details.subjects(index),'megname') && ~isempty(aap.acq_details.subjects(index).megname))
                    tmp=aap.acq_details.subjects(index).megname;
                else
                    numpos=findstr('CBU',aap.acq_details.subjects(index).mriname);
                    tmp=aap.acq_details.subjects(index).mriname(numpos(1):length(aap.acq_details.subjects(index).mriname));
                    tmp=strtok(tmp,' /\\_,.');
                end
                directory= tmp;
            case 2
                directory=sprintf('S%02d',index);
            case 0
                directory=aap.directory_conventions.subject_directory_names{index};
            case 3
                directory=aap.acq_details.subjects(index).mriname;
            otherwise
                aas_log(1,sprintf('Unknown subject directory format (aap.directory_conventions.subject_directory_format=%d',aap.directory_conventions.subject_directory_format));
        end;
end;