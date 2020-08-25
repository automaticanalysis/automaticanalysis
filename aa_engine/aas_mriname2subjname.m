function subjName = aas_mriname2subjname(aap,name)

subjName = sprintf(aap.directory_conventions.subjectoutputformat,name);

subjName=strtok(subjName,' */\\_,.');
