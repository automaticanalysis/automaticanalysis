function [aap passed]=aas_checkpath(aap,pth,nme,allowwildcards)

if (~exist('allowwildcards','var'))
    allowwildcards=false;
end;

if (isempty(pth))
    passed=true;
else
    if (allowwildcards)
        matches=regexp(pth,'[-*?\w/.\_]*','match');
    else
        matches=regexp(pth,'[-\w/.\_]*','match');
    end;
    
    if (length(matches)==1) && (strcmp(matches{1},pth))
        passed=true;
    else
        passed=false;
        aas_log(aap,true,sprintf('Paths in aa can only contain the characters a-z, A-Z, 0-9, _, -, ., /\nYour path %s=''%s'' does not satisfy this.',nme,pth));
    end;
end
