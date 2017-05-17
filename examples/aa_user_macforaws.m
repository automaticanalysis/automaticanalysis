% Automatic analysis
% User master script
% Rhodri Cusack MRC CBU Cambridge Jan 2006-Feb 2010

% This is an example of a branched analysis specified from within a user script
aap=aarecipe('aap_parameters_defaults.xml','aap_tasklist_typical_fmri.xml');

aap=aas_localconfig(aap);

% DEFINE STUDY SPECIFIC PARAMETERS
aap.options.aa_minver=1.0; % will only work on aa version 1.0 or above
 % The study analysis directory 
if (ismac)
 aap.acq_details.root = '/Users/rhodri/testaa';
aap.directory_conventions.rawdatadir='/Users/rhodri/rawmridata';
else
 aap.acq_details.root = '/home/rhodri/testaa';
aap.directory_conventions.rawdatadir='/home/rhodri/rawmridata';
end;

aap.directory_conventions.analysisid='analysis1234';

aap.options.autoidentifystructural=0;
aap.options.autoidentifyfieldmaps=0;

 aap.tasksettings.aamod_slicetiming.sliceorder=[1:36];

% The subjects
aap=aas_addsubject(aap,'Patient_CBU080001_CBU080001/Visit_20080102_092821.203000',[5 6]);

aap.acq_details.subjects(1).structural=3;
aap=aas_addsession(aap,'vstm');
aap=aas_addsession(aap,'esta');

aap.directory_conventions.outputformat='splitbymodule';

aap.acq_details.s3.bucket='awsaa.rhodri.raw';
aap.acq_details.s3.root='processeddata_new';
aap.directory_conventions.rawdatadir='Machine_MRC-CBU_MRC35119_TrioTim';
aap.directory_conventions.seriesoutputformat='Series_%04d';
aap.directory_conventions.remotefilesystem='s3'; % none, s3
aap.options.wheretoprocess='localsingle';   %aws, localsingle, localparallel
aap.directory_conventions.dicomfilter='Image';

aap.options.copystructuraltocentralstore=0;

% SET ANY OTHER PARAMETERS YOU WOULD LIKE TO BE DIFFERENT FROM THE DEFAULTS

% DO PROCESSING

save('~/testaa/myaap.mat','aap')
%[s w]=system('`which aa_doprocessing_parallel` ~/testaa/myaap.mat');

%if (s)
%    aas_log(aap,true,w);
%end;
%aa_doprocessing(aap);

