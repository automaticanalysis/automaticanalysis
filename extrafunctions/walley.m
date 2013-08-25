% Index of the (first) minimum
% Tibor Auer MRC CBU Cambridge 2012-2013

function res = walley(a)
  for i = 1:length(a)
    if a(i) == min(a), break, end
  end
  res = i;
end