function [s, w]=aas_runfslcommand(aap,fslcmd,ENV)

% Setup
fslsetup=deblank(aap.directory_conventions.fslsetup);
if not(isempty(fslsetup))
    if not(fslsetup(end)==';'), fslsetup=[fslsetup ';']; end
    fslsetup = ['source ' fslsetup];
end

if nargin < 3, ENV = {}; end

[s, SPMtool] = aas_cache_get(aap,'spm');
if ~s, aas_log(aap,true,'SPM is not found'); end

% there is a glitch under OS-X wherein Matlab doesn't properly
% inherit $PATH depending on how it was started (e.g., icon double
% click vs command line). As such, runfslcommand can crash because the 
% process spawned by aas_shell can't find FSL. A simple workaround is 
% to define a $PATH entry here in the environment variable passed to 
% aas_shell that explicitly includes "<fsldir>/bin" :

if (aas_ismac)
    ENV = vertcat(ENV,{...
        'PATH', [aap.directory_conventions.fsldir '/bin:$PATH']; ...
        'FSLDIR', aap.directory_conventions.fsldir; ...
        'FSLOUTPUTTYPE', aap.directory_conventions.fsloutputtype; ...
        'MATLABPATH', SPMtool.toolPath;...
        });
else
    ENV = vertcat(ENV,{...
        'FSLDIR', aap.directory_conventions.fsldir; ...
        'FSLOUTPUTTYPE', aap.directory_conventions.fsloutputtype; ...
        'MATLABPATH', SPMtool.toolPath;...
    });
end

switch (aap.directory_conventions.fslshell)
    case 'none'
        [s, w]=aas_shell(fslsetup);
        for e = 1:size(ENV,1)
            setenv(ENV{e,1},ENV{e,2});
        end
        [s, w]=aas_shell(fslcmd);
    case {'csh' 'tcsh'}
        cmd=[aap.directory_conventions.fslshell ' -c ''' fslsetup];
        for e = 1:size(ENV,1)
            cmd = [cmd sprintf('setenv %s %s;',ENV{e,1},ENV{e,2})];
        end
        cmd = [cmd fslcmd ''''];
        [s, w]=aas_shell(cmd);
    case { 'bash' 'zsh' }
        cmd=[aap.directory_conventions.fslshell ' -c ''' fslsetup];
        for e = 1:size(ENV,1)
            cmd = [cmd sprintf('export %s=%s;',ENV{e,1},ENV{e,2})];
        end
        cmd = [cmd fslcmd ''''];
        [s, w]=aas_shell(cmd);
    case 'wsl'
        cmd=['bash' ' -c -i ''' fslsetup];
        fslcmd = convertWinPathToWSL(fslcmd);
        for e = 1:size(ENV,1)
            cmd = [cmd sprintf('export %s=%s;',ENV{e,1},convertWinPathToWSL(ENV{e,2}))];
        end
        cmd = [cmd fslcmd ''''];
        [s, w]=aas_shell(cmd);
end

% Display error if there was one
if (s)
    aas_log(aap,true,sprintf('Error running %s, which was %s',cmd,w));
end