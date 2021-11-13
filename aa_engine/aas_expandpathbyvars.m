% expand environment variables (anything with a $ prefix) in the input path x
% with their current value. This assumes bash-style expansion, where variables are
% referenced $likethis and expressions are evaluated $(like this). We also support
% historical csh-style expansion where expressions are evaluated `like this`.
%
% If x is a cell, struct or multi-row char array we recurse and return
% something of similar structure.
%
% If verbose is true, we print each expanded string to the command window.
%
% Examples:
% aas_expandpathbyvars('/imaging/$USER/aa'); % /imaging/jc01/aa
% aas_expandpathbyvars('$HOME/aa/$HOSTNAME'); % /home/jc01/login24
% aas_expandpathbyvars(aap); % yes, this works
%
% 20180730 J Carlin
%
% xnew = aas_expandpathbyvars(x, [verbose=false])
function x = aas_expandpathbyvars(x, verbose)

if ispc()
    % Not supported yet for windows
    % See issues #288 and #289
    return
end

if ~exist('verbose','var') || isempty(verbose)
    verbose = false;
end

if isstruct(x)
    for thisfn = fieldnames(x)'
        fnstr = thisfn{1};
        for n = 1:numel(x)
            x(n).(fnstr) = aas_expandpathbyvars(x(n).(fnstr), verbose);
        end
    end
    return
end

if iscell(x)
    for cellind = 1:numel(x)
        x{cellind} = aas_expandpathbyvars(x{cellind}, verbose);
    end
    return
end

if size(x,1) > 1
    for rowind = 1:size(x,1)
        x(rowind,:) = aas_expandpathbyvars(x(rowind,:), verbose);
    end
    return
end

% so if we get here it's a single-row array, presumably string
% (nb this also means we ignore other type fields in e.g. struct inputs)
if ischar(x)
    if ~isempty(expandind(x))
        xold = x;
        [err, x] = aas_shell(['echo ' x]);
        % white space and row breaks aren't valid XML so these are errors we assume
        x = deblank(x);
        if isempty(x) || err~=0
            aas_log([],true,sprintf(...
                'failed to expand x: %s',...
                xold));
        end
        if verbose
            aas_log([],false,sprintf(...
                'expanded %s to produce %s',...
                xold, x));
        end
    end
end

function ind = expandind(x)

ind = sort([strfind(x,'$'),strfind(x,'`')]);
