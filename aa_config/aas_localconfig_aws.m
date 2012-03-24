function aap=aas_localconfig_aws(aap)

% Where to process
aap.options.wheretoprocess='aws';   %aws, localsingle, localparallel

% Where to store data 
aap.directory_conventions.remotefilesystem='s3'; % options: none, s3

