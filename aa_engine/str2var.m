function val=str2var(str)
% Can this string be converted to a number? if so than do it.
val = str;
if (numel(str)==0), return; end
digits = '[Inf,NaN,pi,\t,\n,\d,\+,\-,\*,\.,e,i, ,E,I,\[,\],\;,\,]';
s = regexprep(str, digits, ''); % remove all the digits and other allowed characters
if (~all(~isempty(s)))          % if nothing left than this is probably a number
  str(strcmp(str,'\n')) = ';';  % parse data tables into 2D arrays, if any
  try                           % try to convert to a date, like 2007-12-05
    datenum(str);               % if successful than leave it alone
  catch                         % if this is not a date than ...
    num = str2num(str);         % ... try converting to a number
    if(isnumeric(num) && numel(num)>0), val=num; end % if a number than save
  end
end