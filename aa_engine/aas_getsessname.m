function [nme]=aas_getsessname(aap,j)

sessions = aap.acq_details.([aas_getsesstype(aap) 's']);

if ~isempty(sessions(j).name)
    nme = sessions(j).name;
else
    nme='(unknown)';
end;
