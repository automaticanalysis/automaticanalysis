% Automatic Analysis function to retrieve image lists
% function [imagefns]=aas_getfiles_bystream(aap,[i,[j,]]streamname[,inputmodulenumber])
%  stream may be at study, subject or session level depending on number of
%  parameters
%  streamname is name of stream
%
% See also aas_getfiles_bystream, which is intended for EPI images. 
%
% Rhodri Cusack MRC CBU Cambridge, Feb 2010
% Tibor Auer MRC CBU Cambridge, 2012-2013

function [allfiles md5]=aas_getfiles_bystream(aap,varargin)

[inpstreamdesc localroot]=aas_getinputstreamfilename(aap,varargin{:});

if (~exist(inpstreamdesc,'file')) % Try output [TA]
    [inpstreamdesc localroot]=aas_getoutputstreamfilename(aap,varargin{:});
end
if (~exist(inpstreamdesc,'file'))
    aas_log(aap,true,sprintf('Attempting to load stream from file %s, but not found',inpstreamdesc));
end;
fid=fopen(inpstreamdesc,'r');

ind=0;

% There should be an MD5 at the top
lne=fgetl(fid);
if ((length(lne)>3) && strcmp(lne(1:3),'MD5'))
    allfiles=[];
    md5=lne;
else
    allfiles=lne;
    ind=ind+1;
end;


% Now read in the files
while (~feof(fid))
    ind=ind+1;
    allfiles=strvcat(allfiles,fullfile(localroot,fgetl(fid)));
end;
fclose(fid);
end