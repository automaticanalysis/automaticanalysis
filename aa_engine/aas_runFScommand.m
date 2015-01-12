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
setenv('PATH',pth);

setenv('FREESURFER_HOME',aap.directory_conventions.freesurferdir);
setenv('FSF_OUTPUT_FORMAT','nii.gz');

%setenv('MNI_PERL5LIB', fullfile(aap.directory_conventions.freesurferdir, 'mni/lib/perl5/5.8.5'))
%setenv('PERL5LIB', '/opt/freesurfer/mni/lib/perl5/5.8.5')

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

switch (aap.directory_conventions.freesurfershell)
    case {'none','tcsh','bash'}
        cmd=[FSsetup FScmd];
    case 'csh'
        cmd=[FSsetup 'csh -c ' FScmd];
end

disp(cmd)
[s w]=aas_shell(cmd);

% Display error if there was one
if (s)
    aas_log(aap,false,sprintf('Error running %s, which was %s',cmd,w));
end;