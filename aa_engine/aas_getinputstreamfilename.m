% Automatic Analysis function to get input filename for a stream
% Rhodri Cusack Western University May 2012

function [inpstreamdesc localroot]=aas_getinputstreamfilename(aap,varargin)

streamname=varargin{end};
if isstruct(streamname), streamname = streamname.CONTENT; end

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

MAXIT = 1;
if isfield(aap.options,'maximumretry'); MAXIT = aap.options.maximumretry; end
it = 0;

while it < MAXIT
    % Get stream name
    inpstreamdesc=fullfile(localroot,sprintf('stream_%s_inputto_%s.txt',streamname,aap.tasklist.currenttask.name));
    
    % If that doesn't exist...
    if (~exist(inpstreamdesc,'file'))
        % ...look for fully qualified inputs
        fqi_filter=fullfile(localroot,sprintf('stream_*.%s_inputto_%s.txt',streamname,aap.tasklist.currenttask.name));
        fqi=dir(fqi_filter);
        if (numel(fqi)>1)
            aas_log(aap, true, sprintf('Found more than one stream matching filter %s - try fully qualifying stream inputs in this module?',fqi_filter));
        elseif numel(fqi)==1
            inpstreamdesc=fullfile(localroot,fqi(1).name);
            break;
        elseif isempty(fqi)
            it = it + 1;
            pause(0.1);
        end
    else
        break
    end
end
end