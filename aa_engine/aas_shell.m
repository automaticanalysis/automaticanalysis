% aa subroutine - wrapper for shell call
%
% [s,w]=aas_shell(cmd,quiet,stopforerrors)
function [s,w]=aas_shell(cmd,quiet,stopforerrors)

if nargin < 2, quiet=true; end
if nargin < 3, stopforerrors=true; end

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
    elseif ispc()
        % TBD: Color suppression needed on Windows? And how to? For now,
        % just keep empty.
        prefix='';
    else
        prefix='setenv TERM dumb;';
    end

    % cache back
    aas_cache_put([],'shellprefix',prefix);
end

%% Run
log_path = cmd;
if ispc()
    log_path = strrep(log_path,'\','\\'); % If Windows, fix the path for cprintf
end
if ~quiet, aas_log([],false,['Running: ' prefix log_path]); end
[s,w]=system([prefix cmd]);

if strcmp('shell-init: error', w)
    % ensure this is an error code
    s = 1;
end

% ln returns s=1 on "file exists" -- this shouldn't be treated as
% an error when using the -f option (and aas_retrieve_intputs does)

if (s && contains(cmd,'ln')  && contains(cmd,'-f') && contains(w,'File exists'))
    s = 0;
    w = '';      
end


%% Process error if we're in non-quiet mode OR if we want to stop for errors
if s && (~quiet || stopforerrors)
    [junk, wenv]=system('/usr/bin/env');
    aas_log([],false,sprintf('***LINUX ERROR FROM SHELL %s\n***WHILE RUNNING COMMAND\n%s',w,[prefix log_path]));
    aas_log([],stopforerrors,sprintf('***WITH ENVIRONMENT VARIABLES\n%s',wenv));
    aas_log([],false,'***END, CONTINUING');
end

rw = ''; lw = '';
if ~isempty(w)
    l = textscan(w,'%s','delimiter',char(10)); l = l{1};
    % strip off tcsh: errors at start (see http://www.mathworks.com/support/solutions/data/1-18DNL.html?solution=1-18DNL)
    toRem = cell_index(l,'mathworks'); if any(toRem), l(toRem) = []; end
    % select the last shell error
    toShow = find(cell_index(l,'tcsh:') | cell_index(l,'/bin/sh:') | cell_index(l,'/bin/bash:'),1,'last');
    if ~isempty(toShow)
        l = l(toShow:end);
        lw = l{1}; l(1) = [];
    end
    rw = sprintf('%s\n',l{:});
end
w = sprintf('%s\n%s',rw,lw);
if ispc()
    w = strrep(w,'\','\\'); % If Windows, fix the path for logging
end
if ~quiet, aas_log([],false,w); end

