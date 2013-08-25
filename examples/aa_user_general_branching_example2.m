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
 aap.acq_details.root = '/Users/rhodri/testaa_branchexample2';
aap.directory_conventions.rawdatadir='/Users/rhodri/rawmridata';
else
 aap.acq_details.root = '/home/rhodri/testaa_branchexample2';
aap.directory_conventions.rawdatadir='/home/rhodri/rawmridata';
end;

aap.directory_conventions.analysisid='analysis_branch_example3';

aap.options.copystructuraltocentralstore=0;

 aap.tasksettings.aamod_slicetiming.sliceorder=[1:36];

% The subjects
aap=aas_addsubject(aap,'*CBU080001/*',[5 6]);

aap=aas_addsession(aap,'vstm');
aap=aas_addsession(aap,'esta');

aap.directory_conventions.outputformat='splitbymodule';

aap.directory_conventions.remotefilesystem='none'; % none, s3
aap.options.wheretoprocess='localsingle';   %aws, localsingle, localparallel

% Add any extra tasks
extraparameters=[];
extraparameters.aap.tasklist.currenttask.settings.FWHM=14;
extraparameters.aap.directory_conventions.analysisid_suffix='_smoothFWHM14';
[aap]=aas_addtask(aap,'aamod_smooth',[],'aamod_normwrite',extraparameters)

extraparameters=[];
extraparameters.aap.tasklist.currenttask.settings.FWHM=16;
extraparameters.aap.directory_conventions.analysisid_suffix='_smoothFWHM16';
[aap]=aas_addtask(aap,'aamod_smooth',[],'aamod_normwrite',extraparameters)

% SET ANY OTHER PARAMETERS YOU WOULD LIKE TO BE DIFFERENT FROM THE DEFAULTS

% DO PROCESSING

aa_doprocessing(aap);

