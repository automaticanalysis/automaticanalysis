% Automatic Analysis function to retrieve image lists
% function [imagefns]=aas_getfiles_bystream(aap,[i,[j,]]streamname[,inputmodulenumber])
%  stream may be at study, subject or session level depending on number of
%  parameters
%  streamname is name of stream
%
% See a;sp aas_getfiles_bystream, which is intended for EPI images. 
%
% Rhodri Cusack MRC CBU Cambridge, Feb 2010

function [allfiles md5]=aas_getfiles_bystream(aap,varargin)

streamname=varargin{end};

streamname=aas_remapstreamname(aap,streamname,true);

switch (nargin)
    case 2
        localroot=aas_getstudypath(aap);
    case 3
        localroot=aas_getsubjpath(aap,varargin{1});
    case 4
        localroot=aas_getsesspath(aap,varargin{1},varargin{2});
end;

allfiles=[];

% Load in all images by stream name
inpstreamdesc=fullfile(localroot,sprintf('stream_%s_inputto_%s.txt',streamname,aap.tasklist.currenttask.name));

if (~exist(inpstreamdesc,'file'))
    % look for fully qualified inputs
    fqi_filter=fullfile(localroot,sprintf('stream_*.%s_inputto_%s.txt',streamname,aap.tasklist.currenttask.name));
    fqi=dir(fqi_filter);
    if (length(fqi)>1)
        aas_log(aap, true, sprintf('Found more than one stream matching filter %s - try fully qualifying stream inputs in this module?',fqi_filter));
    elseif (length(fqi)==1)
        inpstreamdesc=fullfile(localroot,fqi(1).name);
    end;
end;

if (~exist(inpstreamdesc,'file'))
    aas_log(aap,true,sprintf('Attempting to load stream %s from file %s, but not found',streamname,inpstreamdesc));
end;
fid=fopen(inpstreamdesc,'r');

ind=0;

% There should be an MD5 at the top
lne=fgetl(fid);
if ((length(lne)>3) && strcmp(lne(1:3),'MD5'))
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