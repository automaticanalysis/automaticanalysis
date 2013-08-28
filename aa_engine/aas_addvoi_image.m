% Automatic analysis - add a voi for use by aap_vois_extract
% space of image must be same as that of EPIs

function [aap]=aas_addvoi_image(aap,name,filename) 

newvoi=[];
newvoi.type='image';
newvoi.name=name;
newvoi.filename=filename;

aap.spmanalysis.vois=[aap.spmanalysis.vois {newvoi}];
