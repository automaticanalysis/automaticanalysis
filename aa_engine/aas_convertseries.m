% Automatic Analysis DICOM to NIFTI convertor
% Uses SPM to convert from DICOM to NIFTI-1
%  i=subject number
%  seriesnum=series number
%  outputpath=place to write NIFTI
% Rhodri Cusack MRC CBU, Cambridge 2005

% One of:
%  aas_convertseries_fromstream(aap,inputstream,outputstream[,outputpathsuffix])
%  aas_convertseries_fromstream(aap,i,inputstream,outputstream[,outputpathsuffix])
%  aas_convertseries_fromstream(aap,i,j,inputstream,outputstream[,outputpathsuffix])
%  streamname=inputstreamname
%  i=subject number; j=session number
%  seriesnum=series number
%  outputpath=place to write NIFTI
% Rhodri Cusack MRC CBU, Cambridge 2005
% 
function [aap out dicomheader]=aas_convertseries_fromstream(aap,varargin)
v=varargin;
if (length(v)>3 && ischar(v{end-2}))
    outputpathsuffix=v{end};
    v(end)=[];
else
    outputpathsuffix='';
end;

inputstream=v{end-1};
outputstream=v{end};
v(end-1:end)=[];
    
switch(length(v))
    case 0
        outputpath=aas_getstudypath(aap);
    case 1
        i=varargin{1};
        outputpath=aas_getsubjpath(aap,i);
    case 2
        i=varargin{1};
        j=varargin{2};
        outputpath=aas_getsesspath(aap,i,j);
    otherwise
        aas_log(aap,true,'internal error - wrong number of arguments to aas_convertseries_fromstream.');
end

outputpath_withsuffix=fullfile(outputpath,outputpathsuffix);

aap=aas_makedir(aap,outputpath_withsuffix);

currdir=pwd;
cd(outputpath);

dicomdata=aas_getfiles_bystream(aap,v{:},inputstream);

if (isempty(dicomdata))
    aas_log(aap,1,sprintf('Did not find a dicom series called %s',dicomdirsearchpth));
end;

dicomsearchpth=fullfile(aap.directory_conventions.rawdatadir,aap.acq_details.subjects(i).mriname,dicomdatadir.name,aap.directory_conventions.dicomdatafilter);
dicomdata=dir(dicomsearchpth);


if (length(dicomdata)==0)
    aas_log(aap,1,sprintf('Did not find any dicom data (%s) in %s',aap.directory_conventions.dicomdatafilter,dicomsearchpth));
end;

% Limit number of volumes read in at a time
chunksize_volumes=16;
k=1;

currdir=pwd;
cd (outputpath);
while (k<=size(dicomdata,1))
%    fprintf('***New chunk\n');
    oldAcquisitionNumber=-1;
    thispass_numvolumes=0;
    DICOMHEADERS=[];
    while (k<=size(dicomdata,1))
        DICOMHEADERS=[DICOMHEADERS spm_dicom_headers(fullfile(outputpath,dicomdata(k,:)))];
        if (DICOMHEADERS{end}.AcquisitionNumber~=oldAcquisitionNumber)
            thispass_numvolumes=thispass_numvolumes+1;
            if (thispass_numvolumes>chunksize_volumes)
                DICOMHEADERS=DICOMHEADERS(1:(end-1));
                break;
            end;
        end;
        oldAcquisitionNumber=DICOMHEADERS{end}.AcquisitionNumber;
%        fprintf('Acq %d slice %d\n',oldAcquisitionNumber,k);
        k=k+1;
    end;
    if (~exist('echonumbers','var'))
        DICOMHEADERS_selected=DICOMHEADERS;
    else
        DICOMHEADERS_selected=[];
        for l=1:length(DICOMHEADERS);
            if(any(DICOMHEADERS{l}.EchoNumbers==echonumbers))
                DICOMHEADERS_selected=[DICOMHEADERS_selected DICOMHEADERS(l)];
            end;
        end;
    end;
    conv=spm_dicom_convert(DICOMHEADERS_selected,'all','flat','nii');
    out=[out conv.files];  
end;
out=unique(out); 
cd (currdir);

dicomheader=DICOMHEADERS;

