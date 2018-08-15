% Download the minimal AA demo dataset from the Cusack lab's AWS store - or just
% reference it in aap if it's already available.
%
% aap = aa_downloaddemo(aap)
function aap = aa_downloaddemo(aap)

% where to look for raw DICOM data
sources = strsplit(aap.directory_conventions.rawdatadir,':')';
demoind = cell_index(sources, 'aa_demo');
assert(~any(demoind==0), ...
    'did not find demo directory in aap.directory_conventions.rawdatadir')
assert(numel(demoind)==1,...
    'multiple demo directories in aap.directory_conventions.rawdatadir')
% something like [root]/aa_demo/rawdata
rawdir = sources{demoind};
if exist(rawdir,'dir')
    % we're done
    return
end

% we need to download
rawdirroot = fileparts(fileparts(rawdir));
if ~exist(rawdirroot,'dir')
    success = mkdir(rawdirroot);
    assert(success);
end
aas_log(aap, false, ['downloading demo data to ' rawdirroot]);
% attempt to download the data
outfn = websave(fullfile(rawdirroot,'aa_demo_data.tar.gz'),...
    'http://cusacklab.s3.amazonaws.com/html/downloads/aa_demo_v1.tar.gz');
s = aas_shell(sprintf('tar -xvf %s -C %s',outfn,rawdirroot));
assert(s==0);
% the sub-directory should now exist
assert(exist(sources{demoind},'dir')~=0);
aas_log(aap, false, 'done');
