% Like matlab fileparts, except that '.nii.gz' counts as a single extension
function [pth nme ext]=aas_fileparts(fullpth)

[pth nme ext]=fileparts(deblank(fullpth));

if ~isempty(nme)
    [pth2 nme2 ext2]=fileparts(nme);

    if strcmp(ext2,'.nii') && strcmp(ext,'.gz')
        nme=nme2;
        ext=[ext2 ext];
    end;
end;