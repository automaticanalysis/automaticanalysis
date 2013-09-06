% DICOM locator
%
% 90952,4 -->
% /mridata/cbu/CBU090952_MR09032/20090828_131456/Series_004_CBU_EPI_BOLD_260/1st.dcm
% Tibor Auer MRC CBU Cambridge 2012-2013

function strDcm = mri_finddcm(aap,vol,ser)

strSubj = mri_findvol(aap,vol,1);
strDcmDir = dir(fullfile(strSubj,sprintf(aap.directory_conventions.seriesoutputformat,ser)));
strDcm = dir(fullfile(strSubj,strDcmDir(1).name,aap.directory_conventions.dicomfilter));
strDcm = fullfile(strSubj,strDcmDir(1).name,strDcm(1).name);

