% Automatic Analysis function to retrieve image lists
% function [imagefns]=aas_getfiles_bystream(aap,[i,[j,]]streamname[,inputmodulenumber][,priority])
%  stream may be at study, subject or session level depending on number of
%  parameters
%  streamname is name of stream
%  priority is either input or output
%
% See also aas_getfiles_bystream, which is intended for EPI images. 
%
% Rhodri Cusack MRC CBU Cambridge, Feb 2010
% Tibor Auer MRC CBU Cambridge, 2012-2013

function [allfiles md5]=aas_getfiles_bystream(aap,varargin)

order = {'input' 'output'};
if strcmp(varargin{end},'input') || strcmp(varargin{end},'output')
    if strcmp(varargin{end},'output')
        order = {'output' 'input'};
    end
    varargin(end) = []; 
end

for e = 1:numel(order)
    eval(sprintf('[inpstreamdesc localroot]=aas_get%sstreamfilename(aap,varargin{:});',order{e}));
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
end