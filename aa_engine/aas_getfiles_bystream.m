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

function [allfiles md5 inpstreamdesc]=aas_getfiles_bystream(aap,varargin)

% Indices are all the numeric inputs
reqestedIndices = [varargin{cellfun(@(x) isnumeric(x), varargin)}];

% We need the domain of the stream that we are fetching
if (nargin > 2) && ischar(varargin{1}) % unless the streamname is the only input
    streamDomain = varargin{1}; % nice if it's passed in...
else
    switch numel(reqestedIndices)
        
        case 0
            streamDomain = 'study';
        case 1
            streamDomain = 'subject';
        case 2
            switch aas_getmodality(aap)
                case 'FMRI'
                    streamDomain = 'session';
                    % for backwad compatibilty
                    if ~isempty(strfind(aap.tasklist.currenttask.name,'diffusion'))
                        streamDomain = 'diffusion_session';
                    end
                case 'DWI'
                    streamDomain = 'diffusion_session';
                case {'MTI' 'ASL'}
                    streamDomain = 'special_session';
                case 'MEG'
                    streamDomain = 'meg_session';
            end
            varargin = {streamDomain reqestedIndices varargin{cellfun(@(x) ischar(x), varargin)}};
        otherwise
            aas_log(aap, 1, sprintf('Can''t determine the domain for stream ''%s'', givin these indices: %s. Try using this, aas_getfiles_bystream(aap, ''streamDomain'', [%s], ''%s'')', streamName, strjoin(arrayfun(@(x) sprintf('[%d]',x),reqestedIndices, 'UniformOutput', false)), strjoin(arrayfun(@(x) sprintf('%d',x),reqestedIndices, 'UniformOutput', false)), streamName));
    end
end

% Get the domain indices for this level (e.g., stream indices for this subject)
[junk, domainI] = aas_getN_bydomain(aap, streamDomain, reqestedIndices);

order = {'input' 'output'};
if strcmp(varargin{end},'input') || strcmp(varargin{end},'output')
    if strcmp(varargin{end},'output')
        order = {'output' 'input'};
    end
    varargin(end) = [];
end

% If the requested domain indices exist, we try to get the files, otherwise
% we'll just return empty.
if isempty(reqestedIndices) || ismember(reqestedIndices(end), domainI) % allow for domain = study
    
    for e = 1:numel(order)
        if strcmp('input', order{e})
            [inpstreamdesc localroot]=aas_getinputstreamfilename(aap,varargin{:});
        else
            [inpstreamdesc localroot]=aas_getoutputstreamfilename(aap,varargin{:});
        end
        aas_log(aap,false,sprintf('Load stream from file %s...',inpstreamdesc));
        if exist(inpstreamdesc,'file'), break; end
        aas_log(aap,false,sprintf('\b but not found'));
    end
    if (~exist(inpstreamdesc,'file')) 
        aas_log(aap,true,sprintf('Attempting to load stream from file %s, but not found',inpstreamdesc));
        inpstreamdesc = '';
        allfiles = '';
        md5 = '';
        return
    end;
    
    % to avoid conflict in case of parallel execution make a copy of the file in tmp
    tmpfname = tempname;
    nTries = 1;
    while nTries <= aap.options.maximumretry
        try
            copyfile(inpstreamdesc, tmpfname);
            break;
        catch
            nTries = nTries + 1;
        end
    end
    if ~exist(tmpfname,'file')
        aas_log(aap,true,sprintf('ERROR: Error occurred while copying %s to %s',inpstreamdesc,tmpfname));
    end
    
    retries = 0;
    while retries < 5
        % check that file has been written. There may be a slight delay in
        % the file being available.
        if retries > 0
            copyfile(inpstreamdesc, tmpfname);
            pause(1)
        end
        
        fid=fopen(tmpfname,'r');
        lne=fgetl(fid);
        fclose(fid);
        
        if isnumeric(lne)
            retries = retries + 1;
            disp('MD5 line could not be found: waiting 1 second and retrying')            
            pause(1)
            else
            break
        end
    end
    
    % There should be an MD5 at the top
    ind=0;
    fid=fopen(tmpfname,'r');
    lne=fgetl(fid);
    if ((length(lne)>3) && strcmp(lne(1:3),'MD5'))
        allfiles='';
        md5=lne;
    else
        allfiles=lne;
        ind=ind+1;
        if ~exist('md5','var')
            tn = ['/imaging/dp01/temp/debug_' mfilename '.mat'];
            save(tn)
            error(['md5 did not exist: data dump saved to\n /imaging/dp01/temp/debug_' mfilename '.mat'])
        end
    end
    
    % Now read in the files
    while (~feof(fid))
        ind=ind+1;
        allfiles=strvcat(allfiles,fullfile(localroot,fgetl(fid)));
    end;
    fclose(fid);
    
    delete(tmpfname);
else
    inpstreamdesc = '';
    allfiles = '';
    md5 = '';
end

end