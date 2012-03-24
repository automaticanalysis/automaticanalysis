% Automatic Analysis DMDX to .TXT convertor
% Uses awk script to convert from DMDX to .TXT format
%  i=subject number
%  seriesnum=series number
%  outputpath=place to write .TXT format
% Matt Davi, MRC CBU, Cambridge 2006 **INCOMPLETE**
% (based on Rhodri Cusack MRC CBU, Cambridge 2005 "aas_convertseries.m")

function [aap]=aas_convertdmdx(aap,i,seriesnum,outputpath)

% Get zil file name
zildirsearchpth=fullfile(aap.directory_conventions.zildatadir,aap.acq_details.subjects{i},'Raw_Events');
zildatadir=dir(zildirsearchpth);

% check that directory of zil files exists
if (length(zildatadir)==0)
    aas_log(aap,1,sprintf('Did not find a zil directory called %s',zildirsearchpth));
end;

%% START HERE! %%

dicomsearchpth=fullfile(aap.directory_conventions.rawdatadir,aap.acq_details.subjects{i},dicomdatadir.name,'*.dcm');
dicomdata=dir(dicomsearchpth);
if (length(dicomdata)==0)
    aas_log(aap,1,sprintf('Did not find any dicom data (*.dcm) in %s',dicomsearchpth));
end;

for k=1:length(dicomdata)
    DICOMFN=strvcat(DICOMFN,fullfile(aap.directory_conventions.rawdatadir,aap.acq_details.subjects{i},dicomdatadir.name,dicomdata(k).name));
end;
DICOMHEADERS=spm_dicom_headers(DICOMFN);

% Save dicom headers
save(fullfile(outputpath,'dicom_headers'),'DICOMHEADERS');

% Now convert
currdir=cd;
cd(outputpath);
spm_dicom_convert(DICOMHEADERS);
cd(currdir);
