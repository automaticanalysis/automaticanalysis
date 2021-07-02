
% PR testing - dartel using ds000114

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

aap.options.NIFTI4D = 1;

aap.tasksettings.aamod_segment8_multichan.samp=2;
aap.tasksettings.aamod_segment8_multichan.writenormimg=0;
aap.tasksettings.aamod_dartel_normmni.fwhm=1;


% -------------------------------------------------------------------------
% BIDS
% -------------------------------------------------------------------------

% aap.acq_details.input.combinemultiple = true;

aap = aas_processBIDS(aap, [], {'finger_foot_lips'}, {'sub-01'});


% -------------------------------------------------------------------------
% run
% -------------------------------------------------------------------------

aa_doprocessing(aap);

% if directory_conventions.reportname is undefined, skip reporting

if isfield(aap.directory_conventions,'reportname') && ~isempty(aap.directory_conventions.reportname)
    aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
end


aa_close(aap);


