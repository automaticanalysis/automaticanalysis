function v = getmatlabversion
% as double
v = version;
ind = find(v=='(');
v = sscanf(v(1:ind-1),'%d.');
