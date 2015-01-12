function [s w]=aas_runfslcommand(aap,fslcmd)

% Setting paths now done in aa_doprocessing
fslsetup=deblank(aap.directory_conventions.fslsetup);

if not(isempty(fslsetup)) && not(fslsetup(end)==';')
    fslsetup=[fslsetup ';'];
end;

switch (aap.directory_conventions.fslshell)
    case 'none'
        cmd=[fslsetup fslcmd];
        [s w]=aas_shell(cmd);
    case 'csh'
        cmd=['csh -c "' fslsetup 'setenv FSLOUTPUTTYPE ' aap.directory_conventions.fsloutputtype ';' fslcmd '"'];
        [s w]=aas_shell(cmd);
    case 'bash'
        cmd=['bash -c "' fslsetup 'export FSLOUTPUTTYPE=' aap.directory_conventions.fsloutputtype ';' fslcmd '"'];
        [s w]=aas_shell(cmd);
end;

% Display error if there was one
if (s)
    aas_log(aap,false,sprintf('Error running %s, which was %s',cmd,w));
end;