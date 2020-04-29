
% ds002382 (NAMWords userscript)

% ------------------------------------------------------------------------------------------------------------------------------
% INITIALIZATION
% ------------------------------------------------------------------------------------------------------------------------------

more off;
clear all
aa_ver5;

aap = aarecipe('/Users/peellelab/MATLAB_SCRIPTS/aap_parameters_WUSTL.xml','ds002382.xml');

% ------------------------------------------------------------------------------------------------------------------------------
% FSL path init
% ------------------------------------------------------------------------------------------------------------------------------

FSL_binaryDirectory = '/usr/local/fsl/bin'; 
currentPath = getenv('PATH');
if ~contains(currentPath,FSL_binaryDirectory)
    correctedPath = [ currentPath ':' FSL_binaryDirectory ];
    setenv('PATH', correctedPath);
end

% ------------------------------------------------------------------------------------------------------------------------------
% DIRECTORY AND DEFAULTS
% ------------------------------------------------------------------------------------------------------------------------------

aap.acq_details.root = '/Volumes/DATA01/SCRUB_SUBSET';
aap.directory_conventions.analysisid = 'RESULTS_ds002382';
aap.directory_conventions.rawdatadir = '/Volumes/DATA01/NAMWORDS/BIDSCONVERTED';

% PCT if available

% aap.options.wheretoprocess='matlab_pct';
% aap.directory_conventions.poolprofile = 'local';
% aap.options.aaparallel.numberofworkers = 15;

% ------------------------------------------------------------------------------------------------------------------------------
% STREAM CUSTOMIZATION
% ------------------------------------------------------------------------------------------------------------------------------

% sparse sampling adjustment

 if isfield(aap.tasksettings,'aamod_firstlevel_model')
    for index = 1:numel(aap.tasksettings.aamod_firstlevel_model)
            aap.tasksettings.aamod_firstlevel_model(index).xBF.T0 = 1;
    end
 end

% ------------------------------------------------------------------------------------------------------------------------------
% SECOND LEVEL CUSTOMICATION
% ------------------------------------------------------------------------------------------------------------------------------

YOUNGER = {
'sub-01'
'sub-59'
'sub-43'
'sub-58'
'sub-57'
'sub-55'
'sub-56'
'sub-53'
'sub-51'
'sub-50'
'sub-48'
'sub-46'
'sub-44'
'sub-39'
'sub-40'
'sub-41'
'sub-42'
'sub-36'
'sub-37'
'sub-38'
'sub-28'
'sub-29'
'sub-30'
'sub-31'
'sub-32'
'sub-33'
'sub-34'
'sub-35'
'sub-25'
};

OLDER = {
'sub-60'
'sub-61'
'sub-54'
'sub-49'
'sub-47'
'sub-45'
'sub-02'
'sub-03'
'sub-04'
'sub-05'
'sub-06'
'sub-07'
'sub-08'
'sub-09'
'sub-10'
'sub-11'
'sub-12'
'sub-13'
'sub-14'
'sub-15'
'sub-16'
'sub-17'
'sub-18'
'sub-19'
'sub-20'
'sub-21'
'sub-22'
'sub-23'
'sub-24'
'sub-26'
'sub-27'
'sub-52'
};

aap.tasksettings.aamod_secondlevel_randomise2(2).group_one_subjectIDs = YOUNGER;

aap.tasksettings.aamod_secondlevel_randomise2(3).group_one_subjectIDs = OLDER;

aap.tasksettings.aamod_secondlevel_randomise2(4).group_one_subjectIDs = YOUNGER;
aap.tasksettings.aamod_secondlevel_randomise2(4).group_two_subjectIDs = OLDER;

% ------------------------------------------------------------------------------------------------------------------------------
% BIDS input
% ------------------------------------------------------------------------------------------------------------------------------

aap.acq_details.convertBIDSEventsToUppercase = true;
aap.acq_details.omitNullBIDSEvents = true;

aap = aas_processBIDS(aap);

% ------------------------------------------------------------------------------------------------------------------------------
% contrast specification 
% ------------------------------------------------------------------------------------------------------------------------------

aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts_*', '*', 'uniquebysession', '+1xLISTENWORD', 'LW', 'T');
aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts_*', '*', 'uniquebysession', '+1xLISTENNOISE', 'LN', 'T');
aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts_*', '*', 'uniquebysession', '+1xLISTENWORD|-1xLISTENNOISE', 'LW_G_LN', 'T');
aap = aas_addcontrast(aap, 'aamod_firstlevel_contrasts_*', '*', 'uniquebysession', '+1xREPEATCORRECT|-1xREPEATNOISE|-1xLISTENWORD|+1xLISTENNOISE', 'rwGrn_G_lwGln', 'T');
						   
% ------------------------------------------------------------------------------------------------------------------------------
% RUN
% ------------------------------------------------------------------------------------------------------------------------------

aa_doprocessing(aap);
aa_close(aap);

