%function [s,w]=aas_shell(cmd)
% aa subroutine - wrapper for shell call
% 
function [s,w]=aas_shell(cmd,quiet)

if (~exist('quiet','var'))
    quiet=false;
end;

% Cache this as it is remarkably slow
[hit prefix]=aas_cache_get([],'shellprefix');
if ~hit
    [s,wshell]=system('echo $0');
    pos=find(wshell=='/');
    if ~isempty(pos)
        wshell=wshell(pos(end)+1:end);
    end;
    wshell=strtrim(wshell);
    if (strcmp(wshell,'bash') || strcmp(wshell,'sh'))
        prefix='export TERM=dumb;'; % stops colours
    else
        prefix='setenv TERM dumb;'; % stops colours
    end;
    aas_cache_put([],'shellprefix',prefix);
end;



[s,w]=system([prefix cmd]);
if (~quiet)
    if (strcmp('shell-init: error',w))
        aas_log([],false,sprintf('Likely Linux error %s\n',w));
    end;
end;

if (s && ~quiet)
    [s,wenv]=system('/usr/bin/env');
    aas_log([],false,sprintf('***LINUX ERROR FROM SHELL %s\n***WHILE RUNNING COMMAND\n%s\n***WITH ENVIRONMENT VARIABLES\n%s\nEND, CONTINUING\n',w,[prefix cmd],wenv));
end;
% strip off tcsh: errors at start (see http://www.mathworks.com/support/solutions/data/1-18DNL.html?solution=1-18DNL)
[l r]=strtok(w,10);
while ((length(l)>5) && strcmp(l(1:5),'tcsh:')) || ((length(l)>8) && strcmp(l(1:8),'/bin/sh:'))
    [l r]=strtok(r,10);
end;
w=10;
if (~isempty(l)) w=[l w]; end;
if (~isempty(r)) w=[w r]; end;

