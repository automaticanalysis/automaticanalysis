% Check whether tasksettings have changed since the last execution

function out = aas_checktasksettingconsistency(aap,prev_settings,new_settings)
exclude = {'qsub', 'COMMENT', 'timeadded'};

out = true;

fields1 = fieldnames(prev_settings)';
fields2 = fieldnames(new_settings)';

% different fields
incons = setxor(fields1,fields2);
if ~isempty(incons)
    for f = incons
        if any(strcmp(f{1},exclude))
            continue;
        end
        if isfield(prev_settings,f{1})
            aas_log(aap,false,sprintf('INFO: settings <%s> has been removed',f{1}));
            out = isempty(prev_settings.(f{1}));
        end
        if isfield(new_settings,f{1})
            aas_log(aap,false,sprintf('INFO: settings <%s> has been added',f{1}));
            out = isempty(new_settings.(f{1}));
        end
    end
    if ~out, return; end
end

% different values
for f = fields1
    % exclude
    if any(strcmp(f{1},exclude))
        continue;
    end
    % different class
    if ~strcmp(class(prev_settings.(f{1})),class(new_settings.(f{1})))
        aas_log(aap,false,sprintf('INFO: class of <%s> has changed',f{1}));
        out = false;
        return;
    end
    % different size
    if ~all(size(prev_settings.(f{1}))==size(new_settings.(f{1})))
        aas_log(aap,false,sprintf('INFO: size of <%s> has changed',f{1}));
        out = false;
        return;
    end
    if isstruct(prev_settings.(f{1}))
        for i = 1:numel(prev_settings.(f{1}))
            out = out & aas_checktasksettingconsistency(aap,prev_settings.(f{1})(i),new_settings.(f{1})(i));
        end
    end
    % numeric
    if isnumeric(prev_settings.(f{1})) && ~all(prev_settings.(f{1}) == new_settings.(f{1}))
        aas_log(aap,false,sprintf('INFO: value of <%s> has changed',f{1}));
        out = false;
    end
    % string or cell
    if (ischar(prev_settings.(f{1})) || iscell(prev_settings.(f{1}))) && ~all(strcmp(prev_settings.(f{1}),new_settings.(f{1})))
        aas_log(aap,false,sprintf('INFO: value of <%s> has changed',f{1}));
        out = false;
    end
end

end