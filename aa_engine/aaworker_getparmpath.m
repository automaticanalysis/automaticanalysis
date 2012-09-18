% If third parameter is present and true, don't make path

function [pth]=aaworker_getparmpath(aap,aaworkerid,dontmakedir)
if (strcmp(aap.options.wheretoprocess,'aws'))
    aapth=aas_gettempfilename();
else
    username=spm('GetUser');
    if (strcmp(username,'anonymous'))
        [s w]=system('whoami');
        username=deblank(w);
    end;
    %     if (aas_ismac)
    %         aapth=fullfile('/Users',username,'.aa');
    %     else
    %         aapth=fullfile('/home',username,'.aa');
    %     end;
    
    % Following doesn't work when compiled
    % cd ~
    % homePth = pwd;
    % aapth = fullfile(homePth, '.aa');
    
    % this fails on a mac - weird whitespace that is difficult to remove
    %[s w]=aas_shell('cd ~; pwd');
    %aapth=fullfile(deblank(w));
    if ispc
        aapth = getenv('USERPROFILE');
    else
        aapth = getenv('HOME');
    end
    assert(~isempty(aapth),'failed to find home directory');
end;

if (ischar(aaworkerid))
    pth=fullfile(aapth,sprintf('aaworker%s',aaworkerid));
else
    pth=fullfile(aapth,sprintf('aaworker%d',aaworkerid));
end;

makepath=1;
if (exist('dontmakedir','var'))
    if (dontmakedir)
        makepath=0;
    end;
end;

if (makepath)
    aas_makedir(aap,aapth);
    aas_makedir(aap,pth);
end;
