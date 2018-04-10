function peaks = GetPeaks(fn,opt,thresh,sep,save)
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

v = spm_vol(fn);
[m xyz] = spm_read_vols(v);
% m(m<thresh)=0;
m(isnan(m))=0;

loc = round([xyz' ones(size(xyz,2),1)]*inv(v.mat'));

X = 2:size(m,1)-1; 
Y = 2:size(m,2)-1;
Z = 2:size(m,3)-1;

xyz = zeros(numel(X)*numel(Y)*numel(Z), 3);
count = 0;
for ii = 1:length(Z)
    for jj = 1:length(Y)
        for kk = 1:length(X)
            count = count+1;
            xyz(count,1:3) = [X(kk) Y(jj) Z(ii)];
        end
    end
end


V = sub2ind(size(m),xyz(:,1),xyz(:,2),xyz(:,3));
%%% In Plane 
comp = [m(V-1) m(V+1) m(V-size(m,1)) m(V-size(m,1)+1) m(V-size(m,1)-1) m(V+size(m,1)) m(V+size(m,1)+1) m(V+size(m,1)-1)];
%%% In back of plane.
B = (V-(size(m,1)*size(m,2)));
comp = [comp m(B) m(B-1) m(B+1) m(B-size(m,1)) m(B-size(m,1)+1) m(B-size(m,1)-1) m(B+size(m,1)) m(B+size(m,1)+1) m(B+size(m,1)-1)];
%%% In front of plane.
F = (V+(size(m,1)*size(m,2)));
comp = [comp m(F) m(F-1) m(F+1) m(F-size(m,1)) m(F-size(m,1)+1) m(F-size(m,1)-1) m(F+size(m,1)) m(F+size(m,1)+1) m(F+size(m,1)-1)];

%%% In Plane 
indices = [(V-1) (V+1) (V-size(m,1)) (V-size(m,1)+1) (V-size(m,1)-1) (V+size(m,1)) (V+size(m,1)+1) (V+size(m,1)-1)];
%%% In back of plane.
B = (V-(size(m,1)*size(m,2)));
indices = [indices (B) (B-1) (B+1) (B-size(m,1)) (B-size(m,1)+1) (B-size(m,1)-1) (B+size(m,1)) (B+size(m,1)+1) (B+size(m,1)-1)];
%%% In front of plane.
F = (V+(size(m,1)*size(m,2)));
indices = [indices (F) (F-1) (F+1) (F-size(m,1)) (F-size(m,1)+1) (F-size(m,1)-1) (F+size(m,1)) (F+size(m,1)+1) (F+size(m,1)-1)];



c = m(V);
ind = 1:numel(V);%ind = find(c~=0);
L = zeros(numel(ind), size(comp,2));
for ii = 1:26
    if opt == 1
        L(:,ii) = c(ind)>comp(ind,ii);
    end
    
    if opt == -1
        L(:,ii) = c(ind)<comp(ind,ii);
    end
end

if opt==1
    i1 = find(c>=thresh);
end
if opt==-1
    i1 = find(c<=thresh);
end
i2 = find(sum(L(i1,:),2)==26);

peaks = V(i1(i2));
inds = indices(i1(i2),:);

% orig = peaks;
% keyboard;
% %%
% peaks = orig;

if sep > 0;
    go = 1;
    while go==1;
        
        [x y z] = ind2sub(v.dim,peaks);
        tmp = [x y z ones(numel(x),1)]*v.mat';
        D = triu(DistMatF(tmp(:,1:3),3),1);
        D(D==0)=NaN;
        [r c] = find(D==min(D(:)));
        
        if D(r,c)>=sep
            break
        else
            if m(peaks(r))>m(peaks(c))
                x = x(setdiff(1:numel(x),c));
                y = y(setdiff(1:numel(y),c));
                z = z(setdiff(1:numel(z),c));
                peaks = peaks(setdiff(1:numel(peaks),c));
                inds = inds(setdiff(1:size(inds,1),c),:);
            else
                x = x(setdiff(1:numel(x),r));
                y = y(setdiff(1:numel(y),r));
                z = z(setdiff(1:numel(z),r));
                peaks = peaks(setdiff(1:numel(peaks),r));
                inds = inds(setdiff(1:size(inds,1),r),:);
            end
        end
    end
end

if nargin>4
    V = zeros(v.dim);
    V(inds(:))=1;
    v.fname = [v.fname(1:end-4) '_peaks.nii'];
    v.descrip = [];
    v.dt = [2 0];
    spm_write_vol(v,V);
end

