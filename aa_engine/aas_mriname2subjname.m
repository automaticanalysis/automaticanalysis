function subjName = aas_mriname2subjname(name)

subjName = '';

numpos=strfind(name,'CBU');
if ~isempty(numpos)
    subjName =name(numpos(1):length(name));
    subjName=strtok(subjName,' /\\_,.');
end

if isempty(subjName)
    subjName = name;
end