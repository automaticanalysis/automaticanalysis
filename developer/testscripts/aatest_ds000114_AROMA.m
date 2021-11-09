function aatest_ds000114_AROMA(parameterfile, deleteprevious, wheretoprocess)

% This script tests aamod_AROMA_denoise on a single subject using 
% the ds000114 motor task. The data should be downloaded prior
% to running the script
%
% see aatest_ds00114_TEMPLATE.m in $AAHOME/developer for help on
% writing and using a test script

% this requires AROMA-ICA to be installed and added as an aa toolbox
%
% The easiest way to install is prolly:
%
%   % cd /users/abcd1234/tools
%   % sudo git clone https://github.com/maartenmennes/ICA-AROMA.git
%
%   (assuming that the repo link is still valid)
%
% add the correspodning entry to your parameterset
%   <toolbox desc='Toolbox with implemented interface in extrafunctions/toolboxes' ui='custom'>
%       <name desc='Name corresponding to the name of the interface without the "Class" suffix' ui='text'>aroma</name>
%           <dir ui='dir'>/users/abcd1234/tools/ICA-AROMA</dir>
%   </toolbox>

% -------------------------------------------------------------------------
% init
% -------------------------------------------------------------------------

aap = aa_test_inittest(mfilename('fullpath'), parameterfile, deleteprevious, wheretoprocess);

% -------------------------------------------------------------------------
% analysis options
% -------------------------------------------------------------------------

aap.options.autoidentifystructural_choosefirst = 1;
aap.options.autoidentifystructural_chooselast = 0;

aap.options.NIFTI4D = 1;

aap.acq_details.numdummies = 4;	
aap.acq_details.input.correctEVfordummies = 1;

aap.tasksettings.aamod_segment8.writenormimg = 0;
aap.tasksettings.aamod_segment8.samp = 2;
aap.tasksettings.aamod_smooth.FWHM = 5;

% ------------------------------------------------------------------------------------------------------------------------------
% STREAM CUSTOMIZATION
% ------------------------------------------------------------------------------------------------------------------------------

% explicitly feed coregged epi and structural to AROMA

aap = aas_renamestream(aap,'aamod_AROMA_denoise_00001','epi','aamod_coreg_extended_00001.epi');
aap = aas_renamestream(aap,'aamod_AROMA_denoise_00001','structural','aamod_coreg_extended_00001.structural');

% -------------------------------------------------------------------------
% BIDS
% -------------------------------------------------------------------------

aap.acq_details.input.combinemultiple = true;
aap = aas_processBIDS(aap,[],{'finger_foot_lips'},{'sub-01'});

% -------------------------------------------------------------------------
% modeling - contrast specification
% -------------------------------------------------------------------------

aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','sameforallsessions',[1 0 0],'All:Finger','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','sameforallsessions',[0 1 0],'All:Foot','T');
aap = aas_addcontrast(aap,'aamod_firstlevel_contrasts_00001','*','sameforallsessions',[0 0 1],'All:Lips','T');

% -------------------------------------------------------------------------
% run
% -------------------------------------------------------------------------

aa_doprocessing(aap);

% if directory_conventions.reportname is undefined, skip reporting

if isfield(aap.directory_conventions,'reportname') && ~isempty(aap.directory_conventions.reportname)
    aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
end


aa_close(aap);


