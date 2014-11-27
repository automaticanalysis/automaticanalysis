function subjName = aas_megname2subjname(name)

subjName = '';

numpos=strfind(name,'meg');
if ~isempty(numpos)
    subjName =name(numpos(1):length(name));
    subjName=strtok(subjName,'/');
end

if isempty(subjName)
    subjName = name;
end