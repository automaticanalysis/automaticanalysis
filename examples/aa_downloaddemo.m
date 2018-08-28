% Download dataset from the URL into an aa_demo folder as specified in
% aap.directory_conventions.rawdatadir. It also adds data to
% aap.directory_conventions.rawdatadir;
% If no URL is provided then the small aa_demo from Rhodri Cusack will be used.
%
% ds114_test : https://files.osf.io/v1/resources/9q7dv/providers/osfstorage/57e549f9b83f6901d457d162
%
% aap = aa_downloaddemo(aap,[URL, [rawdir]])
function aap = aa_downloaddemo(aap,URL,rawdir)

DEMODIRBASENAME = 'aa_demo';

DATA = {...
    'https://files.osf.io/v1/resources/umhtq/providers/osfstorage/5b7465680b87a00018c6a76f', 'aa_demo';...
    };

if nargin < 2
    URL = DATA{1,1};
    rawdir = DATA{1,2};
end

% where to look for raw DICOM data
% something like [root]/aa_demo
sources = strsplit(aap.directory_conventions.rawdatadir,':')';
demoind = cell_index(sources, DEMODIRBASENAME);
assert(~any(demoind==0), ...
    ['did not find ''' DEMODIRBASENAME ''' directory in aap.directory_conventions.rawdatadir'])
assert(numel(demoind)==1,...
    ['multiple reference to ''' DEMODIRBASENAME ''' directory in aap.directory_conventions.rawdatadir'])
demodir = sources{demoind};
if exist('rawdir','var') && ~isempty(strfind(demodir,fullfile(DEMODIRBASENAME,rawdir))), demodir = fileparts(demodir); end % remove reference

if ~exist('rawdir','var') || ~exist(fullfile(demodir,rawdir),'dir') % we need to download
    aas_makedir(aap,demodir);
    aas_log(aap, false, ['INFO: downloading demo data to ' demodir]);
    % attempt to download the data
    outfn = [tempname '.tar.gz'];
    assert(aas_shell(sprintf('wget -O %s %s',outfn,URL))==0, ...
        ['could not download dataset from ' URL]);
    [s, w] = aas_shell(sprintf('tar -xvf %s -C %s',outfn,demodir));
    assert(s==0);
    % the sub-directory should now exist
    assert(exist(sources{demoind},'dir')~=0);
    aas_log(aap, false, 'INFO: done');
    % delete the downloaded archive
    delete(outfn);
    rawdir = strsplit(w); rawdir = rawdir{1}(1:end-1);
end
% add dataset to rawdatadir
if ~any(strcmp(sources,fullfile(demodir,rawdir)))
    sources{demoind} = fullfile(demodir,rawdir);
    aap.directory_conventions.rawdatadir = strjoin(sources,':');
end