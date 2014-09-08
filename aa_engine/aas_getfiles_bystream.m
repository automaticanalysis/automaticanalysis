% Automatic Analysis function to retrieve image lists
% function [imagefns]=aas_getfiles_bystream(aap,[i,[j,]]streamname[,inputmodulenumber])
%  stream may be at study, subject or session level depending on number of
%  parameters
%  streamname is name of stream
%
% Updated to return empty for indices that are empty, e.g., for subjects
% that are missing a session.
%
% Rhodri Cusack MRC CBU Cambridge, Feb 2010
% Tibor Auer MRC CBU Cambridge, 2012-2013

function [allfiles md5]=aas_getfiles_bystream(aap,varargin)

% Error check: is the stream an input to this module?
streamName = varargin{end};
if isstruct(streamName), streamName = streamName.CONTENT; end
inputStreams = aap.internal.inputstreamsources{aap.tasklist.currenttask.modulenumber}.stream;
streamI = find(strcmp(streamName, {inputStreams.name}), 1, 'first');
if isempty(streamI), aas_log(aap, 1, sprintf('Module %s doesn''t have stream %s as input!',  aap.tasklist.currenttask.name, streamName)); end

% Indices are all the numeric inputs
reqestedIndices = [varargin{cellfun(@(x) isnumeric(x), varargin)}];

% We need the domain of the stream that we are fetching
if ischar(varargin{1})
    streamDomain = varargin{1}; % nice if it's passed in...
else
    switch numel(reqestedIndices)
        
        case 0
            streamDomain = 'study';
        case 1
            streamDomain = 'subject';
        case 2
            streamDomain = 'session';
        otherwise
            aas_log(aap, 1, sprintf('Can''t determine the domain for stream ''%s'', givin these indices: %s. Try using this, aas_getfiles_bystream(aap, ''streamDomain'', [%s], ''%s'')', streamName, strjoin(arrayfun(@(x) sprintf('[%d]',x),reqestedIndices, 'UniformOutput', false)), strjoin(arrayfun(@(x) sprintf('%d',x),reqestedIndices, 'UniformOutput', false)), streamName));
    end
end


% Get the domain indices for this level (e.g., stream indices for this subject)
[~, domainI] = aas_getN_bydomain(aap, streamDomain, reqestedIndices);

% If the requested domain indices exist, we try to get the files, otherwise
% we'll just return empty.
if ismember(reqestedIndices(end), domainI)
    
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
    
else
    allfiles = [];
    md5 = [];
end

end