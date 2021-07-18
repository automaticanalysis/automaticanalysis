
% developer PR test script
%
% description: BIDS multimodal dataset ds000114 -- diffusion
%

% -------------------------------------------------------------------------
% init
% -------------------------------------------------------------------------

clear all;
aa_ver5;

aap = aarecipe([mfilename('fullpath') '.xml']);

% -------------------------------------------------------------------------
% results and data directory specification
% -------------------------------------------------------------------------

temp = split(mfilename,'_');
aap.directory_conventions.analysisid = [ temp{2} '_' temp{3} ];

fprintf('Saving results in: %s/%s\n', aap.acq_details.root, aap.directory_conventions.analysisid);

aap.directory_conventions.rawdatadir = fullfile(aap.directory_conventions.rawdatadir,temp{2});

% -------------------------------------------------------------------------
% analysis options
% -------------------------------------------------------------------------

aap.acq_details.numdummies = 1;
aap.acq_details.input.combinemultiple = 1;
aap.options.autoidentifystructural_choosefirst = 1;

aap.tasksettings.aamod_diffusion_bet.bet_f_parameter = 0.4;

% -------------------------------------------------------------------------
% BIDS
% -------------------------------------------------------------------------

aap = aas_processBIDS(aap,[],[],{'sub-01'});

% -------------------------------------------------------------------------
% run
% -------------------------------------------------------------------------

aa_doprocessing(aap);

% if directory_conventions.reportname is undefined, skip reporting

if isfield(aap.directory_conventions,'reportname') && ~isempty(aap.directory_conventions.reportname)
    aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
end


aa_close(aap);

