% volunteer locator
%   fp - return fullpath
%
% 11 --> meg11_0011_cc710154/110210
% Tibor Auer MRC CBU Cambridge 2012-2013

function strSubj = meeg_findvol(aap,subjpath,varargin)

argParser = inputParser;
argParser.addParameter('fullpath',false,@(x) islogical(x) || isnumeric(x));
argParser.addParameter('date','',@(x) ischar(x) || isstring(x));
argParser.addParameter('probe',false,@(x) islogical(x) || isnumeric(x));
argParser.parse(varargin{:});

fp = argParser.Results.fullpath;
fdate = argParser.Results.date;
probe = argParser.Results.probe;

strSubj = findvol(aap,'meeg',subjpath,fdate,fp,probe);