function nme=aas_getsubjdesc(aap,i)

nme='';
if ~isempty(aap.acq_details.subjects(i).megname)
    nme=[nme 'MEG:' aap.acq_details.subjects(i).megname ' '];
end;
if ~isempty(aap.acq_details.subjects(i).mriname{1})
    if isnumeric(aap.acq_details.subjects(i).mriname)
        nme=[nme 'MRI:' num2str(aap.acq_details.subjects(i).mriname{1})];    
    else
        nme=[nme 'MRI:' aap.acq_details.subjects(i).mriname{1}];
    end
end;
if isempty(nme)
    nme='(unknown)';
end;

nme=['subject "',aap.acq_details.subjects(i).subjname, '" - ', nme];
