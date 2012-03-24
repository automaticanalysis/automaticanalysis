function [aap]=aas_doprocessing_initialisationmodules(aap)


% THE MODULES LISTED IN AAP.TASKLIST.INITIALISATIONMODULES ARE ALWAYS RUN
for k=1:length(aap.tasklist.initialisation.module)
    % allow full path of module to be provided
[stagepath stagename]=fileparts(aap.tasklist.initialisation.module(k).name);

% retrieve description from module
    description=aap.schema.tasksettings.(stagename).ATTRIBUTE.desc;
    % find out whether this module needs to be executed once per study, subject or session
    domain=aap.schema.tasksettings.(stagename).ATTRIBUTE.domain;

    switch (domain)
        case 'study'
            % now run current stage
            aas_log(aap,0,sprintf('INITIALISATION MODULE %s: %s',stagename,description));
            [aap,resp]=aa_feval(fullfile(stagepath,stagename),aap,'doit');
        case 'subject'
            for i=1:length(aap.acq_details.subjects)
                % now run current stage
                aas_log(aap,0,sprintf('INITIALISATION MODULE %s: %s for %s',stagename,description,aas_getsubjname(aap,i)));
                [aap,resp]=aa_feval(fullfile(stagepath,stagename),aap,'doit',i);
            end;
        case 'session'
            for i=1:length(aap.acq_details.subjects)
                for j=aap.acq_details.selected_sessions
                    aas_log(aap,0,sprintf('INITIALISATION MODULE %s: %s for %s',stagename,description,aas_getsessname(aap,i,j)));
                    [aap,resp]=aa_feval(fullfile(stagepath,stagename),aap,'doit',i,j);
                end;
            end;
        otherwise
            aas_log(aap,1,sprintf('Unknown domain %s associated with stage %s',domain,stagename));
    end;
end;
