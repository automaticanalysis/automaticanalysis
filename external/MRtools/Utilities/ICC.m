function [ICC2 MSS MSE MSF df] = ICC(varargin)
%%% ICC(2,1) from Shrout & Fleiss
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

if numel(varargin)>1
    YY = zeros(numel(varargin{1}),numel(varargin));
    for ii = 1:numel(varargin)
        YY(:,ii) = varargin{ii};
    end
else
    YY = varargin{1};
end

n = size(YY,2);
df1 = size(YY,2)-1;
df2 = size(YY,1)-1;
df3 = numel(YY)-df1-df2-1;

% [df1 df2 df3]
SS1 = 0;
YY2 = mean(YY-repmat(mean(YY,2),1,size(YY,2)));
SS1 = SumOfSquares(YY2(:))*size(YY,1);

SS2 = SumOfSquares(mean(YY,2))*n;
YY2 = (YY-repmat(mean(YY,2),1,size(YY,2))-repmat(mean(YY,1),size(YY,1),1));
SS3 = SumOfSquares(YY2(:));

% [SS1 SS2 SS3]

MSF = SS1/df1;
MSS = SS2/df2;
MSE = SS3/df3;

ICC2 =   (MSS-MSE) / (MSS + (n-1)*MSE + (n*(MSF-MSE)/size(YY,1)));

df = [n size(YY,1)];

