function pvpmod(x,parNames)
% PVPMOD - evaluate parameter/value pairs
% pvpmod(x,parNames) assigns the value x(i+1) to the parameter defined by
% the string x(i) in the calling workspace. This is useful to evaluate
% <varargin> contents in an mfile, e.g. to change default settings of any
% variable initialized before pvpmod(x) is called. 
% Before these assignments the code checks whether the inputs are arranged
% in pairs with the first element of each pair being a char or string. If
% optional input argument parNames (a cell array of char arrays) is
% specified, it also ensures that each of the parameter names in x
% are part of the allowable set given in parNames.


% (c) U. Egert 1998
% Modified by H. Hentschke 2017

% new:
% 0. check inputs
if nargin>1
  validateattributes(parNames,{'cell','string'},{'nonempty'});
end
if ~isempty(x)
  % 0. check inputs
  validateattributes(x,{'cell'},{'nonempty'});
  % 1. check that we're really dealing with pairs in x (i.e. varargin has
  % an even number of elements)
  if rem(numel(x),2)~=0
    error('expecting parameter/value pairs')
  end
  for k = 1:2:numel(x)
    % 2. make sure the first input within a pair is a char or string
    validateattributes(x{k},{'char','string'},{'nonempty'});
    if nargin>1
      % 3. make sure strings match
      validatestring(x{k},parNames);
    end
    % now make assignments
    assignin('caller', x{k}, x{k+1});
  end
end