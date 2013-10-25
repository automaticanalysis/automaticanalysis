function subjName = aas_mriname2subjname(mriname)

numpos=strfind(mriname,'CBU');
if ~isempty(numpos)
    subjName=mriname(numpos(1):length(mriname));
else
    subjName = mriname;
end
subjName=strtok(subjName,' /\\_,.');