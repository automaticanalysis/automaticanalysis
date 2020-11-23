function T = aas_inittoolbox(aap,tbxname)
TBX = aap.directory_conventions.toolbox(strcmp({aap.directory_conventions.toolbox.name},tbxname));
switch numel(TBX)
    case 0
        aas_log(aap,true,['ERROR: setting not found for toolbox ' tbxname]);
    case 1
        if ~exist([tbxname 'Class'],'class')
            aas_log(aap,true,['ERROR: no interface class found for toolbox ' tbxname]);
        end
        constr = str2func([tbxname 'Class']); 
        
        params = {};
        if isfield(TBX,'extraparameters') && ~isempty(TBX.extraparameters)
            for p = fieldnames(TBX.extraparameters)
                val = TBX.extraparameters.(p{1});
                if isempty(val), continue; end
                if ischar(val) && contains(val,':'), val = strsplit(val,':'); end
                params{end+1} = p{1};
                params{end+1} = val;
            end
        end
        T = constr(TBX.dir,params{:});

    otherwise
        aas_log(aap,true,['ERROR: more than one settings found for toolbox ' tbxname]);
end