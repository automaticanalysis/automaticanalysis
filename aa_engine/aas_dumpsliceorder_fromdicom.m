function aas_dumpsliceorder_fromdicom(fn)
H=spm_dicom_headers(fn);
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
fprintf('  slice order %s\n',sliceorder);
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