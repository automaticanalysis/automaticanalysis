% Assert that input pth (with fieldname nme) does not contain forbidden
% characters. Called by aas_validatepaths.
%
% [aap, passed] = aas_checkpath(aap, pth, nme, allowwildcards, allowpathsep)
function [aap, passed] = aas_checkpath(aap, pth, nme, allowwildcards, allowpathsep)

passed = true;

if isempty(pth)
    return
end

% Save original pth input for use in log message
% When used in aas_log messages, escape backward slashes from windows paths.
logsafe_path = strrep(pth, '\', '\\');

allowedchars = '-\w/.\_';
if exist('allowwildcards','var') && ~isempty(allowwildcards) && allowwildcards
    allowedchars = [allowedchars, '*?'];
end

if ispc()
    % On windows:
    % - paths can start with a drive letter followed by a : when it is an absolute path
    % - paths natively have a \ as filesep (though can also have /)
    expression = ['([A-Z]:)?[\\', allowedchars, ']*'];
else
    expression = ['[', allowedchars, ']*'];
end

if exist('allowpathsep','var') && ~isempty(allowpathsep) && allowpathsep
    pths = strsplit(pth, pathsep);
else
    pths = {pth};
end

for i = 1:length(pths)
    checkpath = pths{i};
    matches = regexp(checkpath, expression, 'match');
    if length(matches)~=1 || ~(strcmp(matches{1},checkpath))
        passed = false;
        msg = sprintf('Paths in aa can only contain the characters a-z, A-Z, 0-9, _, -, ., / and sometimes wildcards/colons \nYour path %s=''%s'' does not satisfy this.',nme,logsafe_path);
        aas_log(aap, true, msg);
    end
end

end