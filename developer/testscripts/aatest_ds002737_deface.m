function aatest_ds002737_deface(deleteprevious,wheretoprocess)
% developer PR test script
%
% description: deface structural and T2 using Freesurfer
% dataset: ds002737

% -------------------------------------------------------------------------
% init
% -------------------------------------------------------------------------

aap = aa_test_inittest(mfilename('fullpath'),deleteprevious);

% -------------------------------------------------------------------------
% analysis options
% -------------------------------------------------------------------------

aap.options.wheretoprocess = wheretoprocess;

% -------------------------------------------------------------------------
% BIDS
% -------------------------------------------------------------------------

% anatomy data is in session 3
aap.acq_details.input.selected_subjects = {'sub-01'};
aap.acq_details.input.selected_visits = {'ses-03'};

aap = aas_processBIDS(aap);

% -------------------------------------------------------------------------
% run
% -------------------------------------------------------------------------

aa_doprocessing(aap);

% if directory_conventions.reportname is undefined, skip reporting

if isfield(aap.directory_conventions,'reportname') && ~isempty(aap.directory_conventions.reportname)
    aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
end


aa_close(aap);

