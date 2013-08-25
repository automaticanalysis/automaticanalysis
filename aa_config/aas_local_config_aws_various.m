function aap=aas_localconfig(aap)

% no shortening of subject names
aap.directory_conventions.T1template=fullfile('templates','T1.nii');
aap.directory_conventions.seriesnamingconvention='AWS';
aap.directory_conventions.seriesoutputformat='Series_%04d';
aap.directory_conventions.outputformat='splitbymodule';
aap.options.autoidentifyfieldmaps=0;
aap.options.wheretoprocess='aws';   %aws, localsingle, localparallel
aap.directory_conventions.dicomfilter='Image';

% Where to store data 
aap.directory_conventions.remotefilesystem='s3'; % options: none, s3

