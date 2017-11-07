function ismac = aas_ismac
%AAS_ISMAC 1 if running on a Mac, 0 otherwise
%
% Not all Matlab installations have the ISMAC function, provided
% for compatibility.
%
% See also COMPUTER, ISPC.
%
% $Id$

ismac =  ~isempty(strfind(computer, 'MAC'));
