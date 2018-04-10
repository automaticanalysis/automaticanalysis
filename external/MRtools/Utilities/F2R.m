function [r F p] = F2R(type,X,df1,df2)
% N = 32;
% df1 = 1;
% df2 = N-2;
% alpha = .1;
%%% Written by Aaron Schultz (aschultz@martinos.org)
%%% Copyright (C) 2014,  Aaron P. Schultz
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

switch lower(type)
    case 'r'
        r = X;
        r2 = r.^2;  
        F = (r2./df1)./((1-r2)./df2);
        p = 1-cdf('f',F,df1,df2);
    case 'f'
        F = X;
        r = sqrt(1./((1./(F./df2*df1))+1));
        p = 1-cdf('f',F,df1,df2);
    case 'a'
        alpha = X;
        F = icdf('f',1-alpha,1,df2);
        r = sqrt(1/((1/(F/df2*df1))+1));
        p = 1-cdf('f',F,df1,df2);
    otherwise
        error('type must be either R, F, or A');
end



