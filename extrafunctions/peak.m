% Index of the (first) maximum
% Tibor Auer MRC CBU Cambridge 2012-2013

function res = peak(a)
  for i = 1:length(a)
    if a(i) == max(a), break, end
  end
  res = i;
end