function con = mvpaa_balanceCont( con, mode )
%BALANCE_CON Balance a contrast of 1 0 -1 for GLM
%   con => N-D array containing the contrast (univariate = 1D, multivariate
%   = 2D, usually)
%   mode => manner in which to balance the contrasts...
%           0           - keep 0s as baseline (disrupts interval spacing)
%           1 (default) - keep equal intervals, at cost of baseline

if nargin < 2
    mode = 1;
end

if any(con(:) < 0) && any(con(:) > 0)
    if mode == 0
        pos = (con>0).*con; neg = (con<0).*con;
        ratio = sqrt(abs(nansum(pos(:))/nansum(neg(:))));
        pos = pos/ratio; neg = neg*ratio;
        con = pos + neg;
    elseif mode == 1 % default!
        con = con - nanmean(con(:));
    else
        error('There is no such option!')
    end
    con = con./(nansum(abs(con(:)))/2);
end