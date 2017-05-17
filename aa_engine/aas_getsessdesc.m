function [nme]=aas_getsessdesc(aap,i,j)

nme=[aas_getsubjdesc(aap,i),'; session ',aas_getsessname(aap,j)];
