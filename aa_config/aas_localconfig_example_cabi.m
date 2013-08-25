function aap=aas_localconfig(aap)

% CABI specifics
% no shortening of subject names
aap.directory_conventions.subject_directory_format=3;
aap.directory_conventions.subject_filenames_format=4;
aap.options.autoidentifyfieldmaps=0;
aap.directory_conventions.seriesoutputformat='*_%d';
aap.options.copystructuraltocentralstore=0;
aap.directory_conventions.rawdataafterconversionprefix='f';
aap.directory_conventions.protocol_structural='T1_mprage';
aap.directory_conventions.seriesnamingconvention='CABI';
aap.options.copystructuraltocentralstore=false;
aap.directory_conventions.T1template='/home/bitcjshin/spm5/templates/T1.nii';