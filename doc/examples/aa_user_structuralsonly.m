% Automatic analysis
% User master script
% Rhodri Cusack MRC CBU Cambridge 2005

% RESET ALL PARAMETERS
aap=[];

% ANALYSIS RECIPE
%  This one only converts Siemens structurals 
aap=aarecipe_structuralsonly(aap);

% GET ALL THE PARAMETERS FOR THIS RECIPE
aap=aa_init(aap);

% DEFINE STUDY SPECIFIC PARAMETERS
aap.options.aa_minver=1.0; % will only work on aa version 1.0 or above
  % The study directory 
aap.acq_details.root = '/home/rhodri/pvs/cbu';
  % The subjects
aap.acq_details.subjects={'Chaudhry_Amir_Mr_CBU050011/20051130_135933'};

% SET ANY OTHER PARAMETERS YOU WOULD LIKE TO BE DIFFERENT FROM THE DEFAULTS
aap.options.autoidentifyfieldmaps=0;  % don't bother finding fieldmaps as only doing structurals
 
% DO PROCESSING
aa_doprocessing(aap);


