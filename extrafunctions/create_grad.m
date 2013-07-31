% Create an RGB gradient
% Tibor Auer MRC CBU Cambridge 2012-2013

function cmap = create_grad(c1,c2,step)
cmap(1,:) = c1;
cmap(step,:) = c2;
gr = (c2-c1)/(step-2);
c = c1;
for i = 2:step-1
    c = c + gr;
    cmap(i,:) = c;
end
