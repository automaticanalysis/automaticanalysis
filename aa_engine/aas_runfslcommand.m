function [s w]=aas_runfslcommand(aap,fslcmd)

% Setup
fslsetup=deblank(aap.directory_conventions.fslsetup);
if not(isempty(fslsetup))
    if not(fslsetup(end)==';'), fslsetup=[fslsetup ';']; end
    fslsetup = ['source ' fslsetup];
end


ENV = {...
    'MATLABPATH', path;...
    'FSLOUTPUTTYPE', aap.directory_conventions.fsloutputtype ...
    };

switch (aap.directory_conventions.fslshell)
    case 'none'
        [s w]=aas_shell(fslsetup);
        for e = 1:size(ENV,1)
            setenv(ENV{e,1},ENV{e,2});
        end
        [s w]=aas_shell(fslcmd);
    case {'csh' 'tcsh'}
        cmd=[aap.directory_conventions.fslshell ' -c "' fslsetup];
        for e = 1:size(ENV,1)
            cmd = [cmd sprintf('setenv %s %s;',ENV{e,1},ENV{e,2})];
        end
        cmd = [cmd fslcmd '"'];
        [s w]=aas_shell(cmd);
    case 'bash'
        cmd=[aap.directory_conventions.fslshell ' -c "' fslsetup];
        for e = 1:size(ENV,1)
            cmd = [cmd sprintf('export %s=%s;',ENV{e,1},ENV{e,2})];
        end
        cmd = [cmd fslcmd '"'];
        [s w]=aas_shell(cmd);
end;

% Display error if there was one
if (s)
    aas_log(aap,false,sprintf('Error running %s, which was %s',cmd,w));
end;