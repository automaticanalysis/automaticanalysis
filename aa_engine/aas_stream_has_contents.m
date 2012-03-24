% Automatic Analysis function to retrieve image lists
% function [imagefns]=aas_getfiles_bystream(aap,[i,[j,]]streamname[,inputmodulenumber])
%  stream may be at study, subject or session level depending on number of
%  parameters
%  streamname is name of stream
% Rhodri Cusack MRC CBU Cambridge, Feb 2010

function [stream_has_contents]=aas_stream_has_contents(aap,varargin)

streamname=varargin{end};

streamname=aas_remapstreamname(aap,streamname,true);

% If this stream has source, it will be listed in the dependency maps
stream_has_contents=any(strcmp(streamname,{aap.internal.inputstreamsources{aap.tasklist.currenttask.modulenumber}.stream.name}));
    
end