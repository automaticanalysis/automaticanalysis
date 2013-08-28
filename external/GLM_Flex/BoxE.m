function e = BoxE(cv)
%%% Compute Box's Epsilon on the covariance matrix cv.  e is a metric of
%%% non-sphericity.
%%%
%%% Written by Aaron P. Schultz - aschultz@martinos.org
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

k = size(cv,1);
q  = eye(k)-(1/k);

cvp = q*cv;
ee = eig(cvp);
ee = ee(find(ee>0));
e = (sum(ee)^2) / ((k-1)*sum(ee.^2));
