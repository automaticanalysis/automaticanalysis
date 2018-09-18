% Assert that input pth (with fieldname nme) does not contain forbidden
% characters. Called by aas_validatepaths.
%
% [aap passed]=aas_checkpath(aap,pth,nme,allowwildcards,allowcolon)
function [aap passed]=aas_checkpath(aap,pth,nme,allowwildcards,allowcolon)

if ~exist('allowwildcards','var') || isempty(allowwildcards)
    allowwildcards=false;
end;

if ~exist('allowcolon','var') || isempty(allowcolon)
    allowcolon=false;
end;

if (isempty(pth))
    passed=true;
else
    if (allowwildcards)
        if allowcolon
            matches=regexp(pth,'[-*?\w/.\_]*','match');
        else
            matches=regexp(pth,'[-*?\w/.\_:]*','match');
        end;
    else
        if allowcolon
            matches=regexp(pth,'[-\w/.\_:]*','match');
        else
            matches=regexp(pth,'[-\w/.\_:]*','match');
        end;
    end;
    
    if (length(matches)==1) && (strcmp(matches{1},pth))
        passed=true;
    else
        passed=false;
        aas_log(aap,true,sprintf('Paths in aa can only contain the characters a-z, A-Z, 0-9, _, -, ., / and sometimes wildcards/colons \nYour path %s=''%s'' does not satisfy this.',nme,pth));
    end;
end
