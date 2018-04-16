function ndx = sub2ind_aps(siz,loc)
%==============================================================================
%%% This is a rewrite of the built-in matlab function that reworks how
%%% inputs are specified. Point locations are specified as a matrix rather
%%% than as a series of vectors. Comes in very handy if you don't know 
%%% matrix sizes in advance. Try "help sub2ind" for more details.
%%%
%%% Inputs:
%%% siz = matrix size (e.g. size(M))
%%%
%%% loc = a matrix of point locations (same width as number as length of
%%%     siz).
%%%
%%% Output:
%%% ndx = vector indcies for specified matrix subindices.
%%%
%%% Adapted by Aaron Schultz - aschultz@martinos.org

siz = double(siz);
if length(siz)<2
        error('MATLAB:sub2ind:InvalidSize',...
            'Size vector must have at least 2 elements.');
end

v = ones(1,numel(siz));
v(1:length(loc)) = loc;

%Compute linear indices
k = [1 cumprod(siz(1:end-1))];
ndx = 1;
for i = 1:length(v),
    ndx = ndx + (v(i)-1)*k(i);
end
