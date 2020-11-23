% volunteer locator
%   fp - return fullpath
%
% 90952 --> CBU090952_MR09032/20090828_131456
% Tibor Auer MRC CBU Cambridge 2012-2013
%
% CHANGE HISTORY
%
% 05/17 [MSJ] added a try/catch to handle empty directory and OS X .DS_Store weirdness

function strSubj = mri_findvol(aap,subjpath,fp,probe)

if nargin < 4, probe = false; end
if nargin < 3, fp = false; end

strSubj = findvol(aap,'mri',subjpath,'',fp,probe);