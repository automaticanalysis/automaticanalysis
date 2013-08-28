% Automatic analysis
% User master script
% Rhodri Cusack MRC CBU Cambridge Jan 2006-Feb 2010

% This is an example of a branched analysis specified from within a user script
aap=aarecipe('aap_parameters_defaults.xml','aap_tasklist_typical_fmri.xml');

aap=aas_localconfig(aap);

% DEFINE STUDY SPECIFIC PARAMETERS
aap.options.aa_minver=1.0; % will only work on aa version 1.0 or above
 % The study analysis directory 

aap.acq_details.root = '/cn_data/rhodri';

aap.directory_conventions.analysisid='analysis_asn_019';


aap.options.autoidentifyfieldmaps=0;

 aap.tasksettings.aamod_slicetiming.sliceorder=[1:32];

% The subjects
subjlist={14,16,18,20,23,24,25,26,27,28,29};

aap=aas_addsubject(aap,sprintf('Patient_CBU09%04d_.*/.*/',subjlist{i}),[3:11]);

% Manually specify MPRAGE
% aap.options.autoidentifystructural=0;
% aap.acq_details.subjects(1).structural=2;

for sessind=1:9
    aap=aas_addsession(aap,sprintf('session_%d',sessind));
end;

aap.directory_conventions.outputformat='splitbymodule';

aap.acq_details.s3.root='processeddata_new';
aap.directory_conventions.rawdatadir='Machine_MRC-CBU_MRC35119_TrioTim';
aap.directory_conventions.seriesoutputformat='Series_%04d';
aap.directory_conventions.remotefilesystem='s3'; % none, s3
aap.options.wheretoprocess='localsingle';   %aws, localsingle, localparallel
aap.directory_conventions.dicomfilter='Image';

aap.options.copystructuraltocentralstore=0;

% SET ANY OTHER PARAMETERS YOU WOULD LIKE TO BE DIFFERENT FROM THE DEFAULTS

% DO PROCESSING

save('~/myaap.mat','aap')
%[s w]=system('`which aa_doprocessing_parallel` ~/testaa/myaap.mat');

%if (s)
%    aas_log(aap,true,w);
%end;
aa_doprocessing(aap,'rhodri','lorirhodricusacklorirhodricusack','camneuroasn','activestatenetworks')

