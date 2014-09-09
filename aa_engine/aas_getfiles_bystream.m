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

order = {'input' 'output'};
if strcmp(varargin{end},'input') || strcmp(varargin{end},'output')
    if strcmp(varargin{end},'output')
        order = {'output' 'input'};
    end
    varargin(end) = []; 
end

% If the requested domain indices exist, we try to get the files, otherwise
% we'll just return empty.
if ismember(reqestedIndices(end), domainI)

for e = 1:numel(order)
    if strcmp('input', order{e})
		[inpstreamdesc localroot]=aas_getinputstreamfilename(aap,varargin{:});
	else
		[inpstreamdesc localroot]=aas_getoutputstreamfilename(aap,varargin{:});
	end    
	if exist(inpstreamdesc,'file'), break; end
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