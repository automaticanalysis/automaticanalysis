% Returns filename
% Tibor Auer MRC CBU Cambridge 2012-2013
% 2014-02-13: Rhodri Cusack Western, added processing of multi-line inputs

function f = basename(path)

if size(path,1)>1
    f=[];
    for fileind=1:size(path,1)
        [pth nme ext]=fileparts(path(fileind,:));
        f=char(f,nme);
    end;
else
    [p f] = fileparts(path);
end;