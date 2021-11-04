function aatest_ds002737_deface(parameterfile, deleteprevious, wheretoprocess)

% developer PR test script
%
% description: deface structural and T2 using Freesurfer
% dataset: ds002737

% -------------------------------------------------------------------------
% init
% -------------------------------------------------------------------------

aap = aa_test_inittest(mfilename('fullpath'), parameterfile, deleteprevious, wheretoprocess);

% -------------------------------------------------------------------------
% analysis options
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% BIDS
% -------------------------------------------------------------------------

% anatomy data is in session 3
aap = aas_processBIDS(aap,{'ses-03'},[],{'sub-01'});

% -------------------------------------------------------------------------
% run
% -------------------------------------------------------------------------

aa_doprocessing(aap);

% if directory_conventions.reportname is undefined, skip reporting

if isfield(aap.directory_conventions,'reportname') && ~isempty(aap.directory_conventions.reportname)
    aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
end


aa_close(aap);

