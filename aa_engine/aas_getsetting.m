function val = aas_getsetting(aap,settingstring,index)
% Parse setting path
settingpath = textscan(settingstring,'%s','Delimiter','.'); settingpath = settingpath{1};

% Obtain setting
val = aap.tasklist.currenttask.settings;
for f = settingpath'
    if isfield(val,f{1})
        val = val.(f{1});
    else
        aas_log(aap,false,sprintf('WARNING: Setting <%s> is not specified!',settingstring));
        val = [];
    end
end

if ~isempty(val) && (nargin == 3) % index
    if ischar(val)
        val = textscan(val,'%s');
        val = val{1}';
    elseif isnumeric(val) || isstruct(val)
        val = num2cell(val);
    end
    
    try val = val{index};
    catch E
        aas_log(aap,false,sprintf('WARNING (%s): %s requested %s(%d), but', mfilename, E.stack(2).name, settingstring, index));
        aas_log(aap,false,sprintf('WARNING (%s): only %d value(s) are defined in the current settings.', mfilename, numel(val)));
        aas_log(aap,false,sprintf('WARNING (%s): The first value (%0.9g) will be returned instead.', mfilename, val{1}));
        val = val{1};
    end
end
end
