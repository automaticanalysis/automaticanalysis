% aa component used by parallel engine
%  
%  The directory structure is cached on each machine, so when a file
%  is written it doesn't appear for 8-60 seconds on other machines
%  This made flag passing slow. This uses scp to copy a file to the 
%  destination machine so it is available immediately
%
% Rhodri Cusack MRC CBU Aug 2007

function aas_propagateto(destmachine, fn)
    cmd=['scp ' fn ' ' strtrim(destmachine) ':' fn];   
    [s w]=unix(cmd);
    
