function out = ind2sub_aps(siz,ndx)
%IND2SUB Multiple subscripts from linear index.
%   IND2SUB is used to determine the equivalent subscript values
%   corresponding to a given single index into an array.
%
%   [I,J] = IND2SUB(SIZ,IND) returns the arrays I and J containing the
%   equivalent row and column subscripts corresponding to the index
%   matrix IND for a matrix of size SIZ.  
%   For matrices, [I,J] = IND2SUB(SIZE(A),FIND(A>5)) returns the same
%   values as [I,J] = FIND(A>5).
%
%   [I1,I2,I3,...,In] = IND2SUB(SIZ,IND) returns N subscript arrays
%   I1,I2,..,In containing the equivalent N-D array subscripts
%   equivalent to IND for an array of size SIZ.
%
%   Class support for input IND:
%      float: double, single
%
%   See also SUB2IND, FIND.
% 
%   Copyright 1984-2008 The MathWorks, Inc. 
%   $Revision: 1.13.4.4 $  $Date: 2008/05/23 15:34:40 $
%
% 
% adapted by Aaron Schultz, May 5th, 2010 (aschultz@martinos.org)
% The adaptation comes in the form of a single output matrix rather than a
% set of column vector variables for each dimension.
%

nout = max(numel(siz),1);
siz = double(siz);

if length(siz)<=nout,
  siz = [siz ones(1,nout-length(siz))];
else
  siz = [siz(1:nout-1) prod(siz(nout:end))];
end
n = length(siz);
k = [1 cumprod(siz(1:end-1))];
for i = n:-1:1,
  vi = rem(ndx-1, k(i)) + 1;         
  vj = (ndx - vi)/k(i) + 1; 
  out(1:length(vj),i) = vj; 
% %   varargout{i} = vj; 
  ndx = vi;     
end
