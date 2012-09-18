% Automatic Analysis function to get input filename for a stream
% Rhodri Cusack Western University May 2012

function [inpstreamdesc localroot]=aas_getinputstreamfilename(aap,varargin)

streamname=varargin{end};

streamname=aas_remapstreamname(aap,streamname,true);

switch (nargin)
    case 2
        localroot=aas_getstudypath(aap);
    case 3
        localroot=aas_getsubjpath(aap,varargin{1});
    case 4
        if ischar(varargin{1})
            localroot=aas_getpath_bydomain(aap,varargin{1},varargin{2});
        else
            localroot=aas_getsesspath(aap,varargin{1},varargin{2});
        end;
end;

% Get stream name
inpstreamdesc=fullfile(localroot,sprintf('stream_%s_inputto_%s.txt',streamname,aap.tasklist.currenttask.name));

% If that doesn't exist...
if (~exist(inpstreamdesc,'file'))
    % ...look for fully qualified inputs
    fqi_filter=fullfile(localroot,sprintf('stream_*.%s_inputto_%s.txt',streamname,aap.tasklist.currenttask.name));
    fqi=dir(fqi_filter);
    if (length(fqi)>1)
        aas_log(aap, true, sprintf('Found more than one stream matching filter %s - try fully qualifying stream inputs in this module?',fqi_filter));
    elseif (length(fqi)==1)
        inpstreamdesc=fullfile(localroot,fqi(1).name);
    end;
end;
end