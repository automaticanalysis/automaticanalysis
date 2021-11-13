function aatest_ds000114_diffusion(parameterfile, deleteprevious, wheretoprocess)

% developer PR test script
%
% description: BIDS multimodal dataset ds000114 -- diffusion
%

% -------------------------------------------------------------------------
% init
% -------------------------------------------------------------------------

aap = aa_test_inittest(mfilename('fullpath'), parameterfile, deleteprevious, wheretoprocess);

% -------------------------------------------------------------------------
% analysis options
% -------------------------------------------------------------------------

aap.acq_details.numdummies = 1;
aap.options.autoidentifystructural_choosefirst = 1;

aap.tasksettings.aamod_diffusion_bet.bet_f_parameter = 0.8;

% -------------------------------------------------------------------------
% BIDS
% -------------------------------------------------------------------------

% analyising a single session only
aap = aas_processBIDS(aap,{'ses-test'},[],{'sub-01'});

% -------------------------------------------------------------------------
% run
% -------------------------------------------------------------------------

aa_doprocessing(aap);

% if directory_conventions.reportname is undefined, skip reporting

if isfield(aap.directory_conventions,'reportname') && ~isempty(aap.directory_conventions.reportname)
    aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
end


aa_close(aap);
