% In the case where selected sessions is given as a string of session names
%  (e.g., 'attention eyemovements') this parses them into numeric indices
%  of aap.acq_details.sessions
% Rhodri Cusack BMI Western 2016-08-18

function aap=parse_selected_sessions(aap,sessions,varargin)
% Check subselected sessions
selected_sessions=aap.acq_details.selected_sessions;
if ischar(selected_sessions)
    if strcmp(selected_sessions,'*')
        % Wildcard, same as empty
        selected_sessions=1:numel(sessions);
    else
        % Named sessions, parse to get numbers
        sessionnmes = textscan(selected_sessions,'%s'); sessionnmes = sessionnmes{1};
        selected_sessions=[];
        for sessionnme = sessionnmes'
            sessionind = find(strcmp({sessions.name},sessionnme{1}));
            if isempty(sessionind)
                aas_log(aap,true,sprintf('Unknown session %s specified in selected_sessions field of a branch in the tasklist, sessions were %s',sessionnme{1},sprintf('%s ',sessions.name)));
            end;
            selected_sessions=[selected_sessions sessionind];
        end
    end
    
end

if cell_index(varargin,'subject')
    subj = varargin(cell_index(varargin,'subject')+1);
    [junk, subjSess] = aas_getN_bydomain(aap,aas_getsesstype(aap),subj{1});
    selected_sessions = intersect(selected_sessions,subjSess);
end

aap.acq_details.selected_sessions=selected_sessions;

end