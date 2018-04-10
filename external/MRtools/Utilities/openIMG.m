function [M V] = openIMG(fn);
%%% This a simple wrapper for spm and freesurfer functions to read image 
%%% data into the MATLAB workspace.  You will need FS and SPM installations
%%% for this function to work.
%%%
%%% Input:
%%% fn = filename
%%%
%%% Outputs:
%%% M = volume matrix
%%% V = volume headers
%%%
%%% Written by Aaron Schultz, May 5th, 2010 (aschultz@martinos.org)
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

[a b c] = fileparts(strtrim(fn(1,:)));
    
    
if strcmpi('.mgz',c);
    V = MRIread(fn);
    M = V.vol;
    V.vol = [];
    return
end

if strcmpi('.mgh',c);
    V = MRIread(fn);
    M = V.vol;
    V.vol = [];
    return
end

V = spm_vol(fn);
if V(1).dim(2)==1
    %keyboard;
    V = MRIread(strtrim(fn(1,:)));
    M = zeros(size(fn,1),numel(V.vol));
    for ii = 1:size(fn,1)
        V = MRIread(strtrim(fn(ii,:)));
        ss = size(V.vol);
        if numel(ss)<4;
            ss(end+1:4) = 1;
        end
        M(ii,:) = reshape(V.vol, prod(ss(1:3)),ss(4))';
    end
    V.ss = ss;
%     M = V.vol;
    V.vol = [];
else
%     keyboard;
    M = spm_read_vols(V);
end
