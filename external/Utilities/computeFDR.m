function [FDR_alpha FDR_p FDR_t] = computeFDR(fn,df,alpha,pos,mask)
%%% This function is only appropriate for t-maps (e.g. spmT_0001.img)
%%% Compute FDR corrected threshold for a t-map.
%%%
%%% Inputs:
%%% fn = a filename for a t-statistic image.
%%%
%%% df = the degrees of freedom for the test.  If you don't know this
%%% value, it can usually be found in the header.descrip field of the t
%%% image.
%%%
%%% alpha = the desired FDR corrected alpha value.
%%%
%%% pos = a flag to look at either positive or negative values.  If 1 then
%%% positive values will be evaluated if 0 then negative values will be
%%% evaluated
%%%
%%% Outputs:
%%% FDR_alpha = the alpha value that you set
%%%
%%% FDR_p = the corresponding uncorrected p-value threshold for the 
%%% specified FDR alpha value.
%%%
%%% FDR_t = the corresponding uncorrect t-value threshold for the specified
%%% FDR alpha value.
%%%
%%%
%%% Aaron Schultz 05-07-2010 aschultz@martinos.org
%%% Copyright (C) 2011,  Aaron P. Schultz
%%%
%%% Supported in part by the NIH funded Harvard Aging Brain Study (P01AG036694) and NIH R01-AG027435 
%%%
%%% This program is free software: you can redistribute it and/or modify
%%% it under the terms of the GNU General Public License as published by
%%% the Free Software Foundation, either version 3 of the License, or
%%% any later version.
%%% 
%%% This program is distributed in the hope that it will be useful,
%%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%% GNU General Public License for more details.

if nargin==5
    files = char({char(fn) char(mask)});
    v = spm_vol(files);
    m = spm_read_vols(v);
    
    mask = m(:,:,:,end);
    ind = find(mask>0);
    
    mm = [];
    for ii = 1:size(m,4)-1
        tmp = m(:,:,:,ii);
        mm = [mm; tmp(ind)];
    end
else
  if isnumeric(fn)
    mm = fn;
  else
    files = char(fn);
    v = spm_vol(files);
    m = spm_read_vols(v);
    mm = m(:);
  end
end

if numel(df)==2
    pos = 1;
end

if pos==1
    ind = (find(mm > 0));
elseif pos == -1
    ind = (find(mm < 0));
else
    ind = find(abs(mm) > 0);
end

t = abs(mm(ind));
if numel(df)==2
   pp=(1-cdf('f',t,df(1),df(2))); 
else
   pp = (1-cdf('t',t,df));
end
p = sort(pp)';

rh = ((1:numel(p))./numel(p)).*alpha;
ch = p<=rh;
in = (find(ch==1));

% figure(100); plot(sort(t), [p' rh'], 'linewidth',2); shg
figure(100); plot([p' rh'], 'linewidth',2); shg

if p(1)>rh(1)
    in = [];
end

c = 0;
while isempty(in)
    c = c+1;
    if c==1
        disp(['WARNING! No values found for alpha of ' num2str(alpha) '.  Increasing alpha to find the lowest possible FDR correction.']);
    end
    alpha = alpha+.01;
    rh = ((1:numel(p))./numel(p)).*alpha;
    ch = p<=rh;
    in = (find(ch==1));
    figure(100); plot([p' rh'], 'linewidth',2); shg
end

if p(1)>rh(1)
    FDR_alpha = [];
    FDR_p = [];
    FDR_t = [];
    disp('WARNING! No FDR correction possible.');
else
    FDR_alpha = alpha;
    FDR_p = p(in(end));
    t = sort(t,'descend');
    FDR_t = t(in(end));
end
    
