function strCell = strvcat2cell(vCatStr, strtokSwitch)
    if nargin < 2
        strtokSwitch = 1;
    end
    
    if ~ischar(vCatStr)
        error('Input is not a character array')
    end
    if size(vCatStr,1)==1
        if strtokSwitch == 1
            strCell = {};
            while ~isempty(vCatStr)
               [fn, vCatStr] = strtok(vCatStr);
               strCell = [strCell fn];
            end
        else
            strCell = {vCatStr};
        end
    else
        strCell = {};
        for o = 1:size(vCatStr,1)
            strCell = [strCell deblank(vCatStr(o,:))];
        end
    end
end