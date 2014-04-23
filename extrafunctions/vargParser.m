function vargs = vargParser(vargsin, vdefaults, valid)
% vargs = vargParser(vargsin)
%__________________________________________________________________________
%
% Returned Parameters
%
%   vars: a structure with fields named as in your vargsin
%
% Input Parameters
%
%   vargsin   : a matlab varargin cell array, arguments consist of a specifier
%               (name) and a value 
%
%               i.e. {'arg1', 63, 'arg2', 'red'}
%
%   vdefaults : a cell array of the defaults arguments, consists of a specifier
%               (name), a value, and an array of valid values (if a numeric
%               argument) or a cell array of valid strings (if a string
%               argument).
% 
%               i.e. {'arg1', 63, [60:69], 'arg2', 'red', {'red', 'blue'}
%
% Examples
%   vdefaults = {'dostruct', 1, [0 1], 'doEPIs', 1, [0 1], 'method', 'move', {'move' 'copy'}};        
%   vargs = vargParser(varargin, vdefaults);
%__________________________________________________________________________
%   cwild 11/07/2007
% 

    if mod(length(vargsin),2)
        error('Incorrect number of fields in user supplied arguments, must be a multiple of 2 (Specifier + Value). See ''help'' for details.');
    elseif mod(length(vdefaults),3)
        error('Incorrect number of fields in defaults, must be a multiple of 3 (Specifier + Value + Valid Values). See ''help'' for details.');
    end
    
    vargs = struct;   %returning structure
    vdefs = struct;   %temporary storage for valid/default checking
    defaultArgNames = {};
    
    % Load the supplied defaults
    for i = 1 : length(vdefaults)/3
        dName = vdefaults{3*i-2};
        dVal = vdefaults{3*i-1};
        dValid = vdefaults{3*i};
        if ~ischar(dName)
            error(sprintf('Defaults: Name of argument #%d should be a string', i));
        elseif isa(dVal, 'char') && ~isa(dValid, 'cell')
            error(sprintf('Defaults: Argument #%d ''%s'' appears to be a string argument, and valid values should be specified as a cell array of strings', i, dName));
        elseif isa(dVal, 'double') && ~isa(dValid, 'double')
            error(sprintf('Defaults: Argument #%d ''%s'' appears to be a numeric argument, and valid values should be specified as a 1xn matrix of numbers', i, dName));
        end
        vargs.(dName) = dVal;
        vdefs.(dName) = dValid;
        defaultArgNames = [defaultArgNames {dName}];
    end
    
    % Change any values that the user specified
    for i = 1 : length(vargsin)/2
        argName = vargsin{2*i-1};
        argVal = vargsin{2*i};
        
        % Is this argument specified correctly?
        if ~ischar(argName)
            error('User argument name (cell array element %d) should be a string', i);
            
        % Does the specified argument exist in the defaults?    
        elseif ~ismember(argName, defaultArgNames)
            error('Invalid argument, ''%s''', argName);     
        end
        vargs.(argName) = argVal;
    end
    
    % Now check against valid values
    for i = 1 : length(defaultArgNames)
        argName = defaultArgNames{i};
        if length(vdefs.(argName))
            if isa(vargs.(argName), 'char')
                if ~ismember(vargs.(argName), vdefs.(argName))
                    validVals = '';
                    for j = 1: length(vdefs.(argName))
                        validVals = [validVals ',' vdefs.(argName){j}];
                    end
                    validVals = substr(validVals, 1);
                    error(sprintf('Invalid value ''%s'' supplied for ''%s''.\n Valid values are: <%s>', vargs.(argName), argName, validVals));
                end
            elseif isa(vargs.(argName), 'double')
                if ~ismember(vargs.(argName), vdefs.(argName))
                    validVals = '';
                    for j = 1: length(vdefs.(argName))
                        validVals = [validVals ',' num2str(vdefs.(argName)(j))];
                    end
                    validVals = substr(validVals, 1);
                    error(sprintf('Invalid value ''%d'' supplied for ''%s''.\n Valid values are: <%s>', vargs.(argName), argName, validVals));
                end
            end
        end
    end
    
end