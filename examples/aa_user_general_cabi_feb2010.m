% Automatic analysis
% User master script
% Rhodri Cusack MRC CBU Cambridge Jan 2006

% This script does all processing stages, without skull stripping and using
% segment instead of normalize

% ANALYSIS RECIPE
%  This one does everything, without skull stripping, and with segment not
%  normalize 
aap=aarecipe('aap_parameters_defaults.xml','aap_tasklist.xml');

% Any local configuration
aap=aas_localconfig(aap);


% DEFINE STUDY SPECIFIC PARAMETERS
aap.options.aa_minver=3.01; % will only work on aa version 3.01 or above
 % The study analysis directory 
aap.acq_details.root = '/home/jshin/testaa'; 

% Where to find raw data
aap.directory_conventions.rawdatadir='/home/jshin/s_rawdata/cabicore/Rordenbasic';

% Slice order
aap.tasksettings.aamod_slicetiming.sliceorder=[1:2:36 2:2:36]; % 1:2:36 2:2:36 = interleaved ascending


% The subjects
aap=aas_addsubject(aap,'pilotJaemin_20100105',[5 6 7 8]);

aap=aas_addsession(aap,'run1');
aap=aas_addsession(aap,'run2');
aap=aas_addsession(aap,'run3');
aap=aas_addsession(aap,'run4');


% SET ANY OTHER PARAMETERS YOU WOULD LIKE TO BE DIFFERENT FROM THE DEFAULTS


% DO PROCESSING
aa_doprocessing(aap);


