% Check whether tasksettings have changed since the last execution

function out = aas_checktasksettingconsistency(aap,prev_settings,new_settings,varargin)
exclude = {'qsub', 'COMMENT', 'timeadded'};

argParse = inputParser;
argParse.addParameter('settingsRoot','',@ischar);
argParse.addParameter('doAll',false,@(x) islogical(x) || isnumeric(x));
argParse.parse(varargin{:});

if ~isempty(argParse.Results.settingsRoot), settingsroot = [argParse.Results.settingsRoot '.'];
else, settingsroot = ''; end
doall = argParse.Results.doAll;

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
            aas_log(aap,false,sprintf('INFO: settings <%s%s> has been removed',settingsroot,f{1}));
            out = false;
        end
        if isfield(new_settings,f{1})
            aas_log(aap,false,sprintf('INFO: settings <%s%s> has been added',settingsroot,f{1}));
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
        aas_log(aap,false,sprintf('INFO: class of <%s%s> has changed',settingsroot,f{1}));
        out = false;
        if ~doall, return; end
    end
    % different size
    if ~all(size(prev_settings.(f{1}))==size(new_settings.(f{1})))
        aas_log(aap,false,sprintf('INFO: size of <%s%s> has changed',settingsroot,f{1}));
        out = false;
        if ~doall, return; end
    end
    % struct
    if isstruct(prev_settings.(f{1}))
        if ~isstruct(new_settings.(f{1}))
            aas_log(aap,false,sprintf('INFO: <%s%s> is not a struct any more',settingsroot,f{1}));
            out = false;
        else
            for i = 1:numel(prev_settings.(f{1}))
                out = out & aas_checktasksettingconsistency(aap,prev_settings.(f{1})(i),new_settings.(f{1})(i),'settingsRoot',[settingsroot f{1}],'doAll',doall);
            end
        end
    end
    % numeric
    if isnumeric(prev_settings.(f{1})) && ~all(prev_settings.(f{1}) == new_settings.(f{1}),'all')
        aas_log(aap,false,sprintf('INFO: value of <%s%s> has changed',settingsroot,f{1}));
        out = false;
    end
    % string or cell
    if (ischar(prev_settings.(f{1})) || iscell(prev_settings.(f{1}))) && ...
            (numel(prev_settings.(f{1})) ~= numel(new_settings.(f{1})) || ~isequal(prev_settings.(f{1}),new_settings.(f{1})))
        aas_log(aap,false,sprintf('INFO: value of <%s%s> has changed',settingsroot,f{1}));
        out = false;
    end
end

end