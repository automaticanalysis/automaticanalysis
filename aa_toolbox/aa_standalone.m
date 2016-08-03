% Automatic analysis - Entry point for standalone version. 
%
% FORMAT function aa_standalone(fname_tasklist, fname_aa [, 'mridatadir', <dir>][, 'megdatadir', <dir>][, 'connection', <dir>][, 'anadir', <dir>][, 'subj', <subject indices>])
%   Compulsory arguments:
%       - fname_tasklist: aa tasklist (XML file)
%       - fname_aa: pipeline customisation (MATLAB script file editing aap structure)
%   Optional arguments:
%       - 'mridatadir', <dir>: Directory to find raw MRI data (<dir> can be a colon separated list)
%       - 'megdatadir', <dir>: Directory to find raw MEG data
%       - 'connection', <dir>: Directory of a processed pipeline to connect to. If <dir> is a colon separated list (with two elements), then
%            first element is the directory of the the pipeline to connect to
%            second elements is the name of the latest stage to connect to
%       - 'anadir', <dir>: Directory to put analysis to
%       - 'subj', <subject indices>: Index / indices of subject(s) to be included (N.B.: in the order they have been added)

function aa_standalone(fname_tasklist, fname_aa, varargin)

% Load tasklist and customisation
aap=aarecipe('aap_parameters_defaults_CBSU.xml',fname_tasklist);
run(fname_aa)

% Optiona arguments
args = vargParser(varargin);
if isfield(args,'mridatadir'), aap.directory_conventions.rawdatadir = args.mridatadir; end % data
if isfield(args,'megdatadir'), aap.directory_conventions.rawmegdatadir = args.megdatadir; end % data
if isfield(args,'anadir') % output
    aap.acq_details.root = spm_file(args.anadir,'path'); 
    aap.directory_conventions.analysisid = spm_file(args.anadir,'basename'); 
end
if isfield(args,'subj'), aap.acq_details.subjects = aap.acq_details.subjects(args.subj); end % subject selection

if isfield(args,'connection') % connecting
	con = textscan(args.connection,'%s','delimiter',':'); con = con{1};
	if numel(con) == 1, con{2} = ''; % no maxstage specified
    
	aap=aas_doprocessing_initialisationmodules(aap);
    aap.directory_conventions.allowremotecache = 0;
    remotePipes = struct('host',           '', ...
        'directory',   con{1}, ...
        'allowcache',  0, ...
        'maxstagetag', con{2}, ...
        'checkMD5',    1);
    aap=aas_connectAApipelines(aap,remotePipes);
end

%% DO ANALYSIS
aa_doprocessing(aap);
% aa_report(fullfile(aas_getstudypath(aap),aap.directory_conventions.analysisid));