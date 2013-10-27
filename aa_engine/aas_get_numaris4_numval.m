% Taken from spm_dicom_convert
function val = aas_get_numaris4_numval(str,name)
val1 = aas_get_numaris4_val(str,name);
val  = zeros(size(val1,1),1);
for k = 1:size(val1,1)
    val(k)=str2num(val1(k,:));
end;
return;