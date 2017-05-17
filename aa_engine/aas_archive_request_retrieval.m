function [aap s w]=aas_archive_request_retrieval(aap,streampath)
[aap s w]=aas_runpython(aap,sprintf('import aa_archive_queuetool; c=aa_archive_queuetool.aa_archive_retrieve(); c.queueitem("%s");',streampath));
