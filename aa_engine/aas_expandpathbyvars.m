% expand environment variables (anything with a $ prefix) in the input path x
% with their current value.
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
% x = aas_expandpathbyvars(x, [verbose=false])
function x = aas_expandpathbyvars(x, verbose)

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
    while ~isempty(strfind(x, '$'))
        expandind = strfind(x, '$');
        thisind = expandind(1);
        % need this inside loop since x is changing inside the loop
        separators = strfind(x, filesep);
        nextsep = separators(find(separators > thisind, 1, 'first'));
        if isempty(nextsep)
            nextsep = numel(x)+1;
        end
        shellvariable = x(thisind:nextsep-1);
        % need deblank to get rid of possible line break
        [err,res] = aas_shell(...
            ['echo ' shellvariable]);
        res = deblank(res);
        if isempty(res)
            aas_log([],true,sprintf(...
                'failed to expand shell variable in x: %s',...
                shellvariable));
        end
        x = fullfile(x(1:(thisind-1)),res,x(nextsep:end));
        if verbose
            aas_log([],false,sprintf(...
                'expanded %s to produce %s',...
                shellvariable, x));
        end
    end
end
