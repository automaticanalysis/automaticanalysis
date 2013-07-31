function [s w]=aas_runfslcommand(aap,fslcmd)

setenv('FSLDIR',aap.directory_conventions.fsldir);
setenv('FSLOUTPUTTYPE', aap.directory_conventions.fsloutputtype)

% Add to path
[s pth]=system('echo $PATH');
if s
    pth=getenv('PATH');
else
    pth=deblank(pth);
end;
% Check whether ${FSLDIR}/bin is already in there
fslbin=fullfile(aap.directory_conventions.fsldir,'bin');
% Add colons to beginning and end of path in case fslbin is at beginning or
% end and not bracketed by them
sf=strfind([':' pth ':'],[':' fslbin ':']);
if isempty(sf)
    combinedpath=[pth ':' fslbin];
    setenv('PATH',combinedpath);
end;


fslsetup=deblank(aap.directory_conventions.fslsetup);

if not(isempty(fslsetup)) && not(fslsetup(end)==';')
    fslsetup=[fslsetup ';'];
end;
    
switch (aap.directory_conventions.fslshell)
    case 'none'
        cmd=[fslsetup fslcmd];
        [s w]=aas_shell(cmd);
    case 'csh'
        cmd=['csh -c "' fslsetup  fslcmd '"'];
        [s w]=aas_shell(cmd);
    case 'bash'
        cmd=['bash -c "' fslsetup  fslcmd '"'];
        [s w]=aas_shell(cmd);
end;

% Display error if there was one
if (s)
    aas_log(aap,false,sprintf('Error running %s, which was %s',cmd,w));
end;