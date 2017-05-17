function [S0 DT] = dti_nlfit(vdata, bvals, bvecs)
% Nonlienear fit of the DTI tensor on a single voxel timeseries
% FORMAT [S0 DT] = dti_nlfit(vdata, b, bv)
% vdata - single voxel DWI timeseries
% bvals - [1xN] b-values coresponding to vdata with N volumes
% bvecs - [3xN] DW directions coresponding to vdata with N volumes
% S0    - sphericity index
% DT    - [3x3] diffusion tensor (non-diagonalised)
%_______________________________________________________________________
%
% Copyright (C) 2014 MRC Congition and Brain Sciences Unit
%
% Correia MM, Carpenter TA, Williams GB (2009), Looking for the optimal DTI acquisition scheme given a maximum scan time: are more b-values a waste of time?, MRI, 27(2):163-75. 
%
% Marta Correia and Tibor Auer
% $Id: dti_nlfit.m 2014-08-15 16:00:00Z ta02 $


%-This is merely the help file for the compiled routine
error('dti_nlfit.c not compiled')

