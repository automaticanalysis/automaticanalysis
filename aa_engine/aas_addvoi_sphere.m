% Automatic analysis - add a voi for use by aap_vois_extract
% centres and radius are in mm
% function [aap]=aas_addvoi_sphere(aap,name,centres,radius)

function [aap]=aas_addvoi_sphere(aap,name,centres,radius) 

newvoi=[];
newvoi.type='sphere';
newvoi.name=name;
newvoi.centres=centres;
newvoi.radius=radius;

aap.spmanalysis.vois=[aap.spmanalysis.vois {newvoi}];
