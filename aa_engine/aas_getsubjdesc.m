function nme=aas_getsubjdesc(aap,i)

nme='';
for d = 1:numel(aap.acq_details.subjects(i).meegname)
    if ~isempty(aap.acq_details.subjects(i).meegname{d})
        if isnumeric(aap.acq_details.subjects(i).meegname{d})
            nme=[nme '; MEEG:' num2str(aap.acq_details.subjects(i).meegname{d})];
        else
            nme=[nme '; MEEG:' aap.acq_details.subjects(i).meegname{d}];
        end
    end
end
for d = 1:numel(aap.acq_details.subjects(i).mriname)
    if ~isempty(aap.acq_details.subjects(i).mriname{d})
        if isnumeric(aap.acq_details.subjects(i).mriname{d})
            nme=[nme '; MRI:' num2str(aap.acq_details.subjects(i).mriname{d})];
        else
            nme=[nme '; MRI:' aap.acq_details.subjects(i).mriname{d}];
        end
    end
end
if isempty(nme)
    nme='; (unknown)';
end;

nme=['subject "',aap.acq_details.subjects(i).subjname, '" -', nme(2:end)];
