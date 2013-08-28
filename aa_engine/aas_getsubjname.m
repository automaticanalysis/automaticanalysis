function [nme]=aas_getsubjname(aap,i)

nme='';
if (length(aap.acq_details.subjects(i).megname)>0)
    nme=[nme 'MEG:' aap.acq_details.subjects(i).megname ' '];
end;
if (length(aap.acq_details.subjects(i).mriname)>0)
    nme=[nme 'MRI:' aap.acq_details.subjects(i).mriname];
end;
if isempty(nme)
    nme='(unknown)';
end;

nme=['subject ',nme];
