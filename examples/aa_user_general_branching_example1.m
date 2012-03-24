% Automatic analysis
% User master script
% Rhodri Cusack MRC CBU Cambridge Jan 2006-Feb 2010

% ANALYSIS RECIPE
% An example of a branch defined in the tasklist
aap=aarecipe('aap_parameters_defaults.xml','aap_tasklist_branching_example.xml');

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

aap.directory_conventions.analysisid='analysis_branching_example1';

aap.options.copystructuraltocentralstore=0;

aap.tasksettings.aamod_slicetiming(1).sliceorder=[1:36];
aap.tasksettings.aamod_slicetiming(2).sliceorder=[1:36];

% The subjects
aap=aas_addsubject(aap,'*CBU080001/*',[5 6]);

aap=aas_addsession(aap,'vstm');
aap=aas_addsession(aap,'esta');

aap.directory_conventions.outputformat='splitbymodule';

aap.directory_conventions.remotefilesystem='none'; % none, s3
aap.options.wheretoprocess='localsingle';   %aws, localsingle, localparallel


% SET ANY OTHER PARAMETERS YOU WOULD LIKE TO BE DIFFERENT FROM THE DEFAULTS

% DO PROCESSING

aa_doprocessing(aap);


