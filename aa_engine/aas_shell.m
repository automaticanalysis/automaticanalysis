%function [s,w]=aas_shell(cmd)
% aa subroutine - wrapper for shell call
%
function [s,w]=aas_shell(cmd,quiet)

if nargin < 2, quiet=false; end

%% Prepare
% Cache this as it is remarkably slow
[hit, prefix]=aas_cache_get([],'shellprefix');
if ~hit
    % determine active shell
    [junk,wshell]=system('echo $0');
    pos=find(wshell=='/');
    if ~isempty(pos)
        wshell=wshell(pos(end)+1:end);
    end
    wshell=strtrim(wshell);
    
    % stops colours
    if (strcmp(wshell,'bash') || strcmp(wshell,'sh'))
        prefix='export TERM=dumb;';
    else
        prefix='setenv TERM dumb;';
    end
    
    % cache back
    aas_cache_put([],'shellprefix',prefix);
end

%% Run
[s,w]=system([prefix cmd]);

%% Process error
if ~quiet
    if (strcmp('shell-init: error',w))
        aas_log([],false,sprintf('Likely Linux error %s\n',w));
    end
    if s
        [s,wenv]=system('/usr/bin/env');
        aas_log([],false,sprintf('***LINUX ERROR FROM SHELL %s\n***WHILE RUNNING COMMAND\n%s\n***WITH ENVIRONMENT VARIABLES\n%s\nEND, CONTINUING\n',w,[prefix cmd],wenv));
    end
end

rw = ''; lw = '';
if ~isempty(w)
    l = textscan(w,'%s','delimiter',char(10)); l = l{1};
    % strip off tcsh: errors at start (see http://www.mathworks.com/support/solutions/data/1-18DNL.html?solution=1-18DNL)
    toRem = cell_index(l,'mathworks'); if any(toRem), l(toRem) = []; end
    % select the last shell error
    toShow = find(cell_index(l,'tcsh:') | cell_index(l,'/bin/sh:') | cell_index(l,'/bin/bash:'),1,'last');
    l = l(toShow:end);
    if numel(l) > 1, rw = l{2}; end
    if ~isempty(l), lw = l{1}; end
end
w = sprintf('%s\n%s',rw,lw);

