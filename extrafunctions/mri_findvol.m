% CBU volunteer locator
% 90952 --> CBU090952_MR09032/20090828_131456
% Tibor Auer MRC CBU Cambridge 2012-2013

function strSubj = mri_findvol(vol)

params = xml_read('aap_parameters_defaults.xml');
DATA = params.directory_conventions.rawdatadir.CONTENT;
strSubj = dir(fullfile(DATA,sprintf('CBU%06d*',vol)));
if isempty(strSubj)
    error('Subject CBU%06d* not found',vol);
end
strSubjDir = dir(fullfile(DATA,strSubj.name));
strSubj = fullfile(strSubj.name,strSubjDir(3).name);