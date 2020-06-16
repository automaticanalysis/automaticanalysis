function [s, w]=aas_runFScommand(aap,FScmd,quiet)
if nargin < 3, quiet = false; end

% AVG
% Check whether ${FSLDIR}/bin is already in there
pth=getenv('PATH');
FSbin = fullfile(aap.directory_conventions.freesurferdir,'bin');
mniFSbin = fullfile(aap.directory_conventions.freesurferdir,'mni','bin');
fastFSbin = fullfile(aap.directory_conventions.freesurferdir,'fsfast','bin');
% Add colons to beginning and end of path in case FSbin is at beginning or
% end and not bracketed by them
for path2check = {FSbin mniFSbin fastFSbin}
    sf=strfind([':' pth ':'],[':' path2check{1} ':']);
    if (isempty(sf))
        pth = [pth ':' path2check{1}];
    end
end

[s, SPMtool] = aas_cache_get(aap,'spm');
if ~s, aas_log(aap,true,'SPM is not found'); end
ENV = {...
    'PATH',pth;...
    'SUBJECTS_DIR',getenv('SUBJECTS_DIR');...
    'FREESURFER_HOME',aap.directory_conventions.freesurferdir;...
    'FSF_OUTPUT_FORMAT','nii.gz';...
    'MATLABPATH',SPMtool.toolPath;...
    };

FSsetup='';
for setupscript = {deblank(aap.directory_conventions.freesurfersetup) deblank(aap.directory_conventions.freesurferenvironment)}
    if ~isempty(setupscript{1})
        if setupscript{1}(end)~=';', setupscript{1}=[setupscript{1} ';']; end
        FSsetup=[FSsetup 'source ' setupscript{1}];
    end
end

switch aap.directory_conventions.freesurfershell
    case 'none'
        for e = 1:size(ENV,1)
            setenv(ENV{e,1},ENV{e,2});
        end
        cmd = [FSsetup FScmd];
    case 'bash'
        cmd=[aap.directory_conventions.freesurfershell ' -c '''];
        for e = 1:size(ENV,1)
            cmd = [cmd sprintf('export %s=%s;',ENV{e,1},ENV{e,2})];
        end
        cmd = [cmd FSsetup FScmd ''''];
    case {'csh', 'tcsh'}
        cmd=[aap.directory_conventions.freesurfershell ' -c '''];
        for e = 1:size(ENV,1)
            cmd = [cmd sprintf('setenv %s %s;',ENV{e,1},ENV{e,2})];
        end
        cmd = [cmd FSsetup FScmd ''''];
end

aas_log(aap,false,cmd)

[s, w]=aas_shell(cmd,quiet);
if s > 0
    aas_log(aap,true,'freesurfer error');
end

