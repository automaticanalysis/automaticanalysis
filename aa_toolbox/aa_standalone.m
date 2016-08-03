% Automatic analysis - Entry point for standalone version. 
%
% FORMAT function aa_standalone(fname_tasklist, fname_aa [, 'datadir', <dir>][, 'connectdir', <dir>][, 'anadir', <dir>][, 'subj', <subject indices>])
%   - aap: aap structure with parameters and tasklist

function aa_standalone(fname_tasklist, fname_aa, varargin)

% Load tasklist and customisation
aap=aarecipe('aap_parameters_defaults_CBSU.xml',fname_tasklist);
run(fname_aa)

% Optiona arguments
args = vargParser(varargin);
if nargin > 2, aap.directory_conventions.rawdatadir = varargin{1}; end % data
aap.directory_conventions.analysisid = 'release_01'; 
if nargin > 3, aap.acq_details.root = varargin{2}; end % data
if nargin > 2, aap.directory_conventions.rawdatadir = varargin{1}; end % data

aap.acq_details.root = '/imaging/ta02/aa';

aap=aas_doprocessing_initialisationmodules(aap);
aap.directory_conventions.allowremotecache = 0;
remotePipes = struct('host',           '', ...
    'directory',      fullfile(aap.acq_details.root,'Structural'), ...
    'allowcache',     0, ...
    'maxstagetag',   'aamod_dartel_normmni_00001', ...
    'checkMD5',       1);
aap=aas_connectAApipelines(aap,remotePipes);

%% DO ANALYSIS
aa_doprocessing(aap);
% aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));
