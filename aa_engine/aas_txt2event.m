% Short script assumes that the events are stored in the root directory
% under an events folder + the subjects names
function aap = aas_txt2event(aap, modelName, cnames, indx)

if nargin < 4
    indx = [];
end

for p = 1:length(aap.acq_details.subjects)
    for s = 1:length(aap.acq_details.sessions)
        if isempty(indx)
            indx = strfind(aap.acq_details.subjects(p).mriname, '*');
        end
        
        % Get working directory
        if isempty(indx)
            Wdir = fullfile(aap.acq_details.root, 'events', ...
            aap.acq_details.subjects(p).mriname, ...
            aap.acq_details.sessions(s).name);
        else
            Wdir = fullfile(aap.acq_details.root, 'events', ...
            aap.acq_details.subjects(p).mriname(1:indx(1)-1), ...
            aap.acq_details.sessions(s).name);
        end
        
        % Add each event...
        for c = 1:length(cnames)
            ons = load(fullfile(Wdir,['onset_' cnames{c} '.txt']));
            dur = load(fullfile(Wdir,['duration_' cnames{c} '.txt']));
            
            aap = aas_addevent(aap, modelName, ...
                aap.acq_details.subjects(p).mriname, ...
                aap.acq_details.sessions(s).name, ...
                cnames{c}, ...
                ons, dur);
        end
    end
end