function subjName = aas_mriname2subjname(aap,name)

subjName = name;

% subjprefix
subjprefix = aap.directory_conventions.subjectoutputformat;
subjprefix = subjprefix(1:regexp(subjprefix,'[%*]','once')-1);

numpos=strfind(name,subjprefix);
if ~isempty(numpos)
    subjName =name(numpos(1):length(name));
end

subjName=strtok(subjName,' */\\_,.');
