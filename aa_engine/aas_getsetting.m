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
    catch
        aas_log(aap,false,sprintf('WARNING: Setting <%s> has fewer then %d elements.\nWARNING: First value will be applied!',settingstring,index));
        val = val{1};
    end
end
end