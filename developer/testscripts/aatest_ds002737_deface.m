
% developer PR test script
%
% description: deface structural and T2 using Freesurfer
% dataset: ds002737

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

