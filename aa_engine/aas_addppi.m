% Automatic analysis PPI (psychophysiological interaction)
% Rhodri Cusack MRC CBU Cambridge Jan 2008
% Add ppi for processing by aamod_ppi_model

function [aap]=aas_addppi(aap,ppiname, voiname, contrast)

myppi=[];
myppi.ppiname=ppiname;
myppi.voiname=voiname;
myppi.contrast=contrast;

aap.spmanalysis.ppis=[aap.spmanalysis.ppis {myppi}];

