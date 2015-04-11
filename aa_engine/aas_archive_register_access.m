function aas_archive_register(aap,streampath)
if ~isempty(getenv('AA_ARCHIVE'))
    [aap s w]=aas_runpython(aap,sprintf('import aa_archive_queuetool; c=aa_archive_queuetool.aa_archive_access(); c.queueitem("%s");',streampath));
end;