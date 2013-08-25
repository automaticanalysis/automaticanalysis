% Automatic analysis
% User master script
% Rhodri Cusack MRC CBU Cambridge 2005
% ENSURE SYSTEM IN PATH
aa_ver1_devel

% RESET ALL PARAMETERS
aap=[];

% ANALYSIS RECIPE
aap=aarecipe_notfmri_ver01(aap);

% GET ALL THE PARAMETERS FOR THIS RECIPE
aap=aa_init(aap);

% DEFINE STUDY SPECIFIC PARAMETERS
aap.options.aa_minver=1.0; % will only work on aa version 1.0 or above
  % The study directory 
aap.acq_details.root = '/imaging/rhodri/siemens_veins';
  % Add subjects details
   % no EPI, no series to ignore   
   % convert series 12-19 and put them in the specialseries sub directory 
   % MPRAGE and real time t maps will be picked up and converted as usual
aap=aas_addsubject(aap,'*CBU060809/*',[],[],[12:19]); 

% DO PROCESSING
aa_doprocessing(aap);


