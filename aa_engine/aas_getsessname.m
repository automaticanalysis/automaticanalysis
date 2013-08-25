function [nme]=aas_getsessname(aap,i,j)

nme='';
if (length(aap.acq_details.sessions(j).name)>0)
    nme=aap.acq_details.sessions(j).name;
else
    nme='(unknown)';
end;

nme=[aas_getsubjname(aap,i),'; session ',nme];
