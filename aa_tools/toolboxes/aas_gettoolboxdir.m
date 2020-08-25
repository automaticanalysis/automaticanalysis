function tbdir = aas_gettoolboxdir(aap,tbxname)
TBX = aap.directory_conventions.toolbox(strcmp({aap.directory_conventions.toolbox.name},tbxname));
switch numel(TBX)
    case 0
        aas_log(aap,true,['ERROR: setting not found for toolbox ' tbxname]);
    case 1
        tbdir = TBX.dir;
    otherwise
        aas_log(aap,true,['ERROR: more than one settings found for toolbox ' tbxname]);
end