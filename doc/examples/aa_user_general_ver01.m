% Automatic analysis
% User master script
% Rhodri Cusack MRC CBU Cambridge Jan 2006

% This script does all processing stages, without skull stripping and using
% segment instead of normalize

% RESET ALL PARAMETERS
aap=[];

% ANALYSIS RECIPE
%  This one does everything, without skull stripping, and with segment not normalize 
aap=aarecipe_general_ver01(aap);

% MODIFY STANDARD RECIPE MODULE SELECTION HERE IF YOU'D LIKE
  % do this here

% GET ALL THE PARAMETERS FOR THIS RECIPE
aap=aa_init(aap);

% DEFINE STUDY SPECIFIC PARAMETERS
aap.options.aa_minver=1.0; % will only work on aa version 1.0 or above
  % The study directory 
aap.acq_details.root = '/imaging/rhodri/test_segment/aa';
  % The subjects
aap=aas_addsubject(aap,'*CBU050005/*',[4 10]);

aap.acq_details.sessions={'mystery1','mystery2'};
  % Bruker session numbers for EPIs: these correspond to condition names above
%aap.acq_details.brukersessionnums={[4 10],[5 11],[6 12]};
  % Number of dummy scans at the start of each session
aap.acq_details.numdummies=0;


% SET ANY OTHER PARAMETERS YOU WOULD LIKE TO BE DIFFERENT FROM THE DEFAULTS
aap.options.autoidentifyfieldmaps=0;

% To switch back to normalise not segment
%aap.spmanalysis.usesegmentnotnormalise=0;

% DO PROCESSING
aa_doprocessing(aap);


