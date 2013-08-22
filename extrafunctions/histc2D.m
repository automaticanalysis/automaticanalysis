function [p v] = histc2D(r, l, u, m)
if nargin < 4
  m = (r ~= 0);
end
lr = min(r(m));
ur = max(r(m));
dr = ur-lr;
dp = u-l;
if dr
    cr = double(dp)/double(dr);
else
    cr = (l+u)/2*l;
end
for i = l:u
    v(i-l+1) = (i-l)/cr+lr;
end
pr = masked((r-lr)*cr+l,m,false);    
p0 = masked(r,m,true);
p = round(pr + p0);
p(p>u) = u;

