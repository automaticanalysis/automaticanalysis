% If third parameter is present and true, don't make path

function [pth]=aaworker_getparmpath(aap,aaworkerid,dontmakedir)
if strcmp(aap.options.wheretoprocess,'aws')
    aapth=aas_gettempfilename();
elseif ~isempty(aap.options.aaworkerroot)
    aapth = aap.options.aaworkerroot;
else
    username=spm('GetUser');
    if strcmp(username,'anonymous')
        [s, w]=system('whoami');
        username=deblank(w);
    end;
  
    if ispc
        aapth = getenv('USERPROFILE');
    elseif aas_ismac
        aapth=fullfile('/Users',username);
    else
        aapth = getenv('HOME');
    end
    
    if isempty(aapth)
        aapth=fullfile('/home',username);
    end
    
    assert(~isempty(aapth),'failed to find home directory');
end

aapth = fullfile(aapth,'.aa');
if nargin < 2 % return aaworker root
    pth = aapth;
    return
end
    
if ischar(aaworkerid)
    pth=fullfile(aapth,sprintf('aaworker%s',aaworkerid));
elseif isnumeric(aaworkerid)
    pth=fullfile(aapth,sprintf('aaworker%02d',aaworkerid));
end

if (nargin < 3) || ~dontmakedir
    aas_makedir(aap,aapth);
    aas_makedir(aap,pth);
end;
