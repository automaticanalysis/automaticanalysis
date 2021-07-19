% Automatic Analysis function to retrieve image lists
% function [imagefns]=aas_getfiles_bystream(aap,[domain,incides,]streamname)
%  stream may be at study, subject or session level depending on number of
%  parameters
%  streamname is name of stream
% Rhodri Cusack MRC CBU Cambridge, Feb 2010

function [stream_has_contents]=aas_stream_has_contents(aap,varargin)

streamname=varargin{end};

streamname=aas_remapstreamname(aap,streamname,true);

% If this stream has source, it will be listed in the dependency maps
% Test for both simple and fully specified names (if in the depenedency map)
c = cellfun(@(x) [{x} strsplit(x,'.')], {aap.internal.inputstreamsources{aap.tasklist.currenttask.modulenumber}.stream.name},'UniformOutput',false);

stream_has_contents=any(strcmp(streamname,horzcat(c{:})));
    
if numel(varargin) == 3 % domain specific
    aap.options.verbose = -1; % ignore error
    stream_has_contents = ~isempty(aas_getfiles_bystream(aap,varargin{1},varargin{2},streamname,'input'));
end

end