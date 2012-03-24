% Function to check slice order
% Loads up aap_parameters.mat and goes through all sessions for all
% subjects checking dicom_headers.mat
% field CSASeriesInfo, item MrPhoenixProtocol, ASCCONV block,
% sSliceArray.ucMode parameter
% Rhodri Cusack MRC CBU Cambridge May 2008

function aas_dumpsliceorder(fn)

if (isstruct(fn))
    aap=fn;
else
    if (~exist('fn','var'))
        fn='.';
    end;

    load(fullfile(fn,'aap_parameters.mat'));
end;

for i=1:length(aap.acq_details.subjects)
    fprintf('Subject %s\n',aap.acq_details.subjects(i).mriname);
    for j=1:length(aap.acq_details.sessions)
        sesspth=aas_getsesspath(aap,i,j);
        H=load(fullfile(sesspth,'dicom_headers.mat'));
        H=H.DICOMHEADERS(1);        
        ucMode=getpheonix(H,'sSliceArray.ucMode');
        switch ucMode
            case '0x1'
                sliceorder='ascending sequential';
            case '0x2'
                sliceorder='descending sequential';
            case '0x4'
                sliceorder='ascending interleaved';
            otherwise
                sliceorder=sprintf('unrecognised (%s)',ucMode);
        end;
        fprintf('  Session %s slice order %s\n',aap.acq_details.sessions(j).name,sliceorder);
                
    end;
end;

end

function [fieldvalue]=getpheonix(H,fieldname)
fieldvalue='###FIELDNOTFOUND###';
for i=1:length(H{1}.CSASeriesHeaderInfo)
    if (strcmp(H{1}.CSASeriesHeaderInfo(i).name,'MrPhoenixProtocol')) 
        value=H{1}.CSASeriesHeaderInfo(i).item(1).val;
        break
    end;
end;

pos1=strfind(value,'### ASCCONV BEGIN ###');
value=value(pos1:end);
[junk value]=strtok(value,10);
while(length(value)>0)
    [fld value]=strtok(value,10);
    [nme rem]=strtok(fld);
    [eq rem]=strtok(rem);
    if (strcmp(fieldname,nme))
        fieldvalue=strtrim(rem);
        break;
    end;
end;
end
    