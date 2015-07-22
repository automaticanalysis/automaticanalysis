function val = aas_getsetting(aap,settingstring,index)
% Parse setting path
settingpath = textscan(settingstring,'%s','Delimiter','.'); settingpath = settingpath{1};

% Obtain setting
val = aap.tasklist.currenttask.settings;
for f = settingpath'
    val = val.(f{1});
end

if nargin == 3 % index
    if ischar(val)
        val = textscan(val,'%s');
        val = val{1}';
    elseif isnumeric(val)
        val = num2cell(val);
    end
    
    try val = val{index};
    catch
        aas_log(aap,0,sprintf('WARNING: Setting %s has fewer then %d elements.\nWARNING: First value will be applied!',settingstring,index));
        val = val{1};
    end
end
end