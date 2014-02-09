% Returns filename
% Tibor Auer MRC CBU Cambridge 2012-2013

function f = basename(path)

f = [];
for i = 1:size(path,1)
    [p tmp] = fileparts(path(i,:));
    f = strvcat(f,tmp);
end
