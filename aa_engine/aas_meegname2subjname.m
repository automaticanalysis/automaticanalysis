function subjName = aas_meegname2subjname(aap,name)

subjName = sprintf(aap.directory_conventions.meegsubjectoutputformat,name);

subjName=strtok(subjName,'*/');