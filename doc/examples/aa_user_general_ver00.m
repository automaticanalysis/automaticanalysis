% Automatic analysis
% User master script
% Rhodri Cusack MRC CBU Cambridge 2005

% RESET ALL PARAMETERS
aap=[];

% ANALYSIS RECIPE
%  This one only converts Siemens structurals 
aap=aarecipe_general(aap);

% MODIFY STANDARD RECIPE MODULE SELECTION HERE IF YOU'D LIKE
  % do this here

% GET ALL THE PARAMETERS FOR THIS RECIPE
aap=aa_init(aap);

% DEFINE STUDY SPECIFIC PARAMETERS
aap.options.aa_minver=1.0; % will only work on aa version 1.0 or above
  % The study directory 
aap.acq_details.root = '/home/rhodri/pvs/cbu';
  % Add subjects and session numbers for EPI data
aap=aas_addsubject(aap,'*CBU050011/*',[5 11]);
aap=aas_addsubject(aap,'*CBU050012/*',[4 10]);
 % Condition names for each session, must be same for all subjects
aap.acq_details.sessions={'mystery1','mystery2'};
  % Number of dummy scans at the start of each session
aap.acq_details.numdummies=18;


% SET ANY OTHER PARAMETERS YOU WOULD LIKE TO BE DIFFERENT FROM THE DEFAULTS
 
% DO PROCESSING
aa_doprocessing(aap);


