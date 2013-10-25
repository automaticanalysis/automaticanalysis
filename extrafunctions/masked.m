function mm = masked(a, m, i)
  a = double(a);
  if (nargin>2) && i 
    m = not(m);
  end
  mm = a.*m;
