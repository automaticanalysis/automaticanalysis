function C = aas_str2cell(s, delim)
%AAS_STR2CELL Make cell array from delimited string
%
% C = AAS_STR2CELL(S) creates a cell array C from delimited string S.
% For example:
%
%  S = 'one two three';
%  C = STR2CELL(S);
%  % C = {'one' 'two' 'three'}
%
% C = AAS_STR2CELL(S, DELIM) uses the specified delimiter. If nothing is
% specified a comma and then whitespace are tried (either should work).

% First check that string isn't already a cell array. If so, we are done
% here. Otherwise, proceed.
if iscell(s)
  return
end


% If not specified, probably whitespace delimiter but we will check for
% comma. If comma exists, use that; otherwise, use whitespace. Trying to be
% flexible for people's xml files.
if nargin < 2 || isempty(delim)
  if ~isempty(findstr(',',s))
    delim = ',';
  else
    delim = [9:13 32]; % whitespace, copied from STRTOK
  end
end


% Start with an empty cell array and successively append each token,
C = {};
while ~isempty(s)
  [t,s] = strtok(s, delim);
  C{length(C)+1} = strtrim(t); % STRTRIM prevents extra whitespace from sneaking in
end
