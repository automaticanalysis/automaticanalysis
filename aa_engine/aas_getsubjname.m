function nme=aas_getsubjname(aap,i)

nme='';
if ~isempty(aap.acq_details.subjects(i).megname)
    nme=[nme 'MEG:' aap.acq_details.subjects(i).megname ' '];
end;
if ~isempty(aap.acq_details.subjects(i).mriname)
    if isnumeric(aap.acq_details.subjects(i).mriname)
        nme=[nme 'MRI:' num2str(aap.acq_details.subjects(i).mriname)];    
    else
        nme=[nme 'MRI:' aap.acq_details.subjects(i).mriname];
    end
end;
if isempty(nme)
    nme='(unknown)';
end;

nme=['subject ',nme];
