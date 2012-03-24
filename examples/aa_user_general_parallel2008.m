% Automatic analysis
% User master script
% Rhodri Cusack MRC CBU Cambridge Jan 2006

% This script does all processing stages, without skull stripping and using
% segment instead of normalize

% ANALYSIS RECIPE
%  This one does everything, without skull stripping, and with segment not
%  normalize 
aap=aarecipe('aap_parameters_defaults.xml','aap_tasklist.xml');

% DEFINE STUDY SPECIFIC PARAMETERS
aap.options.aa_minver=1.0; % will only work on aa version 1.0 or above
 % The study analysis directory 
aap.acq_details.root = '/imaging/dm01/aaparallel_test'; 
% The subjects
aap=aas_addsubject(aap,{'meg07_0002','*CBU080461/*'},[3 4]);

aap=aas_addsession(aap,'vstm');
aap=aas_addsession(aap,'esta');


% SET ANY OTHER PARAMETERS YOU WOULD LIKE TO BE DIFFERENT FROM THE DEFAULTS


% DO PROCESSING
aa_doprocessing_parallel(aap,'continue');


