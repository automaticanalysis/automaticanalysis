function ind = contains(str, list)
%%% This is a function to find matches to the regular expression defined in  
%%% str within the elements of the cell array list.
%%%
%%% What is returned is a vector of index matches similar to what you get
%%% with the find function only for cell arrays with strings.
%%%
%%% Use .* for wild card
%%% Trailing Parenthese must be specified as \)
%%% ^ means "Starts With
%%%      Metacharacter   Meaning
%%%    ---------------  --------------------------------
%%%                .    Any character
%%%               []    Any character contained within the brackets
%%%              [^]    Any character not contained within the brackets
%%%               \w    A word character [a-z_A-Z0-9]
%%%               \W    Not a word character [^a-z_A-Z0-9]
%%%               \d    A digit [0-9]
%%%               \D    Not a digit [^0-9]
%%%               \s    Whitespace [ \t\r\n\f\v]
%%%               \S    Not whitespace [^ \t\r\n\f\v]
%%%     Metacharacter   Meaning
%%%    ---------------  --------------------------------
%%%              ()     Group subexpression
%%%               |     Match subexpression before or after the |
%%%               ^     Match expression at the start of string
%%%               $     Match expression at the end of string
%%%              \<     Match expression at the start of a word
%%%              \>     Match expression at the end of a word
%%%     Metacharacter   Meaning
%%%    ---------------  --------------------------------
%%%               *     Match zero or more occurrences
%%%               +     Match one or more occurrences
%%%               ?     Match zero or one occurrence
%%%            {n,m}    Match between n and m occurrences
%%%
%%%
%%% Adapted by Aaron Schultz (aschultz@martinos.org) 
%%% Copyright (C) 2011,  Aaron P. Schultz
%%%
%%% Supported in part by the NIH funded Harvard Aging Brain Study (P01AG036694) and NIH R01-AG027435 
%%%
%%% This program is free software: you can redistribute it and/or modify
%%% it under the terms of the GNU General Public License as published by
%%% the Free Software Foundation, either version 3 of the License, or
%%% any later version.
%%% 
%%% This program is distributed in the hope that it will be useful,
%%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%% GNU General Public License for more details.

% str = regexprep(str,')','\\)');
% str = regexprep(str,'(','\\(');

which = regexp(list,str);
ind = [];
for kk = 1:length(which);
    if ~isempty(which{kk})
        ind(end+1) = kk;
    end
end