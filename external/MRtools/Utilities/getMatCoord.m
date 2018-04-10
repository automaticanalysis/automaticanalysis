function [ml vi] = getMatCoord(h,mni,diam)
%%% Get Matrix coordinates and vector index from mni location and a header
%%% structure.
%%% Set diam = (2*voxdim)+1 to plump out the corrdinates by one voxel.
%%% e.g. for 3x3x3 set diam = 7;
%%%
%%% Written by Aaron Schultz (aschultz@martinos.org)
%%%
%%% Last Updated Dec. 11 2012;
%%%
%%% Copyright (C) 2012,  Aaron P. Schultz
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


ind = [];

if nargin == 3 %&& size(mni,1)==1    
    xyz = zeros(prod(h(1).dim),4);
    c = 0;
    for ii = 1:h(1).dim(3)
        for jj = 1:h(1).dim(2)
            for kk = 1:h(1).dim(1)
                c = c+1;
                xyz(c,:) = [kk jj ii 1];
            end
        end
    end
    
    ws = xyz*h(1).mat';
    ws = ws(:,1:3);
    
    ind = [];
    for ii = 1:size(mni,1)
        %keyboard;
        dist = (sqrt(sum((ws-repmat(mni(ii,:),size(ws,1),1)).^2,2)));
        
        adj = 0;
        %tmp = spm_imatrix(h.mat); adj = mean(abs(tmp(7:9)))/2;
        
        ind = [ind; find(dist<=((diam/2)+adj))];
        %numel(find(dist<(diam/2)))
    end
    ind = unique(ind);    
    
    
    ml = round([ws(ind,:) ones(numel(ind),1)]/h(1).mat');
    ml = ml(:,1:3);
    vi = sub2ind(h(1).dim,ml(:,1),ml(:,2),ml(:,3));
end
   
if nargin==2 || isempty(ind);
    mni = [mni ones(size(mni,1),1)];
    ml = round(mni/(h(1).mat'));
    ml = ml(:,1:3);
    vi = sub2ind(h(1).dim,ml(:,1),ml(:,2),ml(:,3));
end