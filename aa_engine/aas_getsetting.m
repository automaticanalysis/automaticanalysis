function [val, index] = aas_getsetting(aap,settingstring,varargin)
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

if ~isempty(val) && (nargin > 2) % index
    if numel(varargin) == 1
        index = varargin{1};
        if ischar(val)
            val = textscan(val,'%s');
            val = val{1}';
        elseif isnumeric(val) || isstruct(val)
            val = num2cell(val);
        end
    else
        switch varargin{1}
            case 'subject'
                index = [];
                if ~isfield(val,'subject'), aas_log(aap,false,'WARNING: there is no subject-specific setting.');
                else
                    index = find(strcmp({val.subject},aas_getsubjname(aap,varargin{2})));
                end
                if isempty(index)
                    aas_log(aap,false,sprintf('WARNING: Setting <%s> for %s is not specified!',settingstring,aas_getsubjdesc(aap,varargin{2})));
                    val = {[]};
                    index = 1;
                end
                if numel(index) > 1
                    aas_log(aap,false,sprintf('WARNING: More than 1 setting <%s> for %s is specified -> only the first will be returned.',settingstring,aas_getsubjdesc(aap,varargin{2})));
                    index = index(1);
                end
             case 'session'
                index = [];
                if ~isfield(val,'session'), aas_log(aap,false,'WARNING: there is no session-specific setting.');
                else
                    index = find(strcmp({val.subject},aas_getsubjname(aap,varargin{2}(1))) & strcmp({val.session},aas_getsessname(aap,varargin{2}(2))));
                end
                if isempty(index)
                    aas_log(aap,false,sprintf('WARNING: Setting <%s> for %s is not specified!',settingstring,aas_getsessdesc(aap,varargin{2}(1),varargin{2}(2))));
                    val = {[]};
                    index = 1;
                end   
        end
        val = num2cell(val);
    end
    
    try val = val{index};
    catch E
        aas_log(aap,false,sprintf('WARNING (%s): %s requested %s(%d), but', mfilename, E.stack(min(2,numel(E.stack))).name, settingstring, index));
        aas_log(aap,false,sprintf('WARNING (%s): only %d value(s) are defined in the current settings.', mfilename, numel(val)));
        aas_log(aap,false,sprintf('WARNING (%s): The first value (%0.9g) will be returned instead.', mfilename, val{1}));
        val = val{1};
    end
end
end
