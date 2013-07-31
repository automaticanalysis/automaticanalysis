% CBU DICOM locator
% 90952,4 -->
% CBU090952_MR09032/20090828_131456/Series_004_CBU_EPI_BOLD_260/1st.dcm
% Tibor Auer MRC CBU Cambridge 2012-2013

function strDcm = mri_finddcm(vol,ser)
params = xml_read('aap_parameters_defaults.xml');
DATA = params.directory_conventions.rawdatadir.CONTENT;
strSubj = dir(fullfile(DATA,sprintf('CBU%06d*',vol)));
if isempty(strSubj)
    error('Subject CBU%06d* not found',vol);
end
strSubjDir = dir(fullfile(DATA,strSubj.name));
strDcmDir = dir(fullfile(DATA,strSubj.name,strSubjDir(3).name,sprintf('Series_%03d*',ser)));
strDcm = dir(fullfile(DATA,strSubj.name,strSubjDir(3).name,strDcmDir(1).name,'*.dcm'));
strDcm = fullfile(DATA,strSubj.name,strSubjDir(3).name,strDcmDir(1).name,strDcm(1).name);

