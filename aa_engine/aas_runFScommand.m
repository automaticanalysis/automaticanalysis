function [s w]=aas_runFScommand(aap,FScmd)
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
    end;
end

ENV = {'PATH',pth;...
    'SUBJECTS_DIR',getenv('SUBJECTS_DIR');...
    'FREESURFER_HOME',aap.directory_conventions.freesurferdir;...
    'FSF_OUTPUT_FORMAT','nii.gz';...
    'MNI_PERL5LIB', fullfile(aap.directory_conventions.freesurferdir, 'mni/lib/perl5/5.8.5');...
    'PERL5LIB', fullfile(aap.directory_conventions.freesurferdir, 'mni/lib/perl5/5.8.5')};

FSsetup='';
for setupscript = {deblank(aap.directory_conventions.freesurfersetup) deblank(aap.directory_conventions.freesurferenvironment)}
    if ~isempty(setupscript{1})
        if setupscript{1}(end)~=';'
            FSsetup=[FSsetup 'source ' setupscript{1} ';'];
        else
            FSsetup=[FSsetup setupscript{1} ];
        end
    end
end

switch aap.directory_conventions.freesurfershell
    case 'none'
        for e = 1:size(ENV,1)
            setenv(ENV{e,1},ENV{e,2});
        end
        cmd = [FSsetup FScmd];
    case 'bash'
        cmd=['bash -c ". ' FSsetup];
        for e = 1:size(ENV,1)
            cmd = [cmd sprintf('export %s=%s;',ENV{e,1},ENV{e,2})];
        end
        cmd = [cmd FScmd '"'];
    case {'csh', 'tcsh'}
        cmd=[FSsetup aap.directory_conventions.freesurfershell ' -c "'];
        for e = 1:size(ENV,1)
            cmd = [cmd sprintf('setenv %s %s;',ENV{e,1},ENV{e,2})];
        end
        cmd = [cmd FScmd '"'];
end

disp(cmd)

[s, w]=aas_shell(cmd);

% Display error if there was one
if (s)
    aas_log(aap,false,sprintf('Error running %s, which was %s',cmd,w));
end;