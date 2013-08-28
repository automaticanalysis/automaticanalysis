function A = catstruct(varargin)
% CATSTRUCT - concatenate structures
%
%   X = CATSTRUCT(S1,S2,S3,...) concates the structures S1, S2, ... into one
%   structure X.
%
%   A.name = 'Me' ; 
%   B.income = 99999 ; 
%   X = CATSTRUCT(A,B) ->
%     X.name = 'Me' ;
%     X.income = 99999 ;
%
%   CATSTRUCT(S1,S2,'sorted') will sort the fieldnames alphabetically.
%
%   If a fieldname occurs more than once in the argument list, only the last
%   occurence is used, and the fields are alphabetically sorted.
%
%   To sort the fieldnames of a structure A use:
%   A = CATSTRUCT(A,'sorted') ;
%
%   See also CAT, STRUCT, FIELDNAMES, STRUCT2CELL

% for Matlab R13
% version 2.0 (sep 2007)
% (c) Jos van der Geest
% email: jos@jasen.nl

% History
% Created:  2005
% Revisions
%   2.0 (sep 2007) removed bug when dealing with fields containing cell arrays (Thanks to Rene Willemink)



N = nargin ;

error(nargchk(1,Inf,N)) ;

if ~isstruct(varargin{end}),
    if isequal(varargin{end},'sorted'),
        sorted = 1 ;
        N = N-1 ;
        if N < 1,
            A = [] ;
            return
        end
    else
        error('Last argument should be a structure, or the string "sorted".') ;
    end
else
    sorted = 0 ;
end

for ii=1:N,
    X = varargin{ii} ;
    if ~isstruct(X),
        error(['Argument #' num2str(ii) ' is not a structure.']) ;
    end
    FN{ii} = fieldnames(X) ;
    VAL{ii} = struct2cell(X) ;
end

FN = cat(1,FN{:}) ;
VAL = cat(1,VAL{:}) ;
[UFN,ind] = unique(FN) ;

if numel(UFN) ~= numel(FN),
    warning('Duplicate fieldnames found. Last value is used and fields are sorted') ;
    sorted = 1 ;
end

if sorted,
    VAL = VAL(ind) ;
    FN = FN(ind) ;
end

% This deals correctly with cell arrays
A = cell2struct(VAL, FN);



