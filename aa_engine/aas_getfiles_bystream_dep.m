function fname = aas_getfiles_bystream_dep(aap,domain,indices,stream)
dep = dep_read(fullfile(aap.internal.aap_initial.acq_details.root,aap.directory_conventions.analysisid,'aap_prov.trp'));
stage = dep.(aap.tasklist.currenttask.name).(stream);

[name, ind] = strtok_ptrn(stage,'_0');
index = sscanf(ind,'_%d');
stages = {aap.tasklist.main.module.name};
istage = find(strcmp(stages,name));
aap = aas_setcurrenttask(aap,istage(index));

fname = aas_getfiles_bystream_multilevel(aap,domain,indices,stream,'output');