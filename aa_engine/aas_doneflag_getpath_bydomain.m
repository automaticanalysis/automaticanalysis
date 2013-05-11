% aa - get path to done flag
% examples:
%  pth=aas_doneflag_getpath(aap,'study',[],4);  % done flag for stage 4
%  pth=aas_doneflag_getpath(aap,'study',[],'s3',4);   % path on s3 
%  pth=aas_doneflag_getpath(aap,'subject',[1 2],4);      % subj 1 sess 2 path for module 4
%  pth=aas_doneflag_getpath(aap,'subject',2,4);      % subj 2 path for module 4
% If not specified, remote filesystem defaults to aap.directory_conventions.remotefilesystem
% This differs from the behaviour of aas_getstudypath and derivations,
%  which default to local filesystem

function [doneflag doneflagpath stagetag]=aas_doneflag_getpath_bydomain(aap,domain,indices,varargin)

%[hit resp]=aas_cache_get(aap,mfilename,domain,indices,varargin{:});
hit=false;
if (hit)
    doneflag=resp{1};
    doneflagpath=resp{2};
    stagetag=resp{3};
else
    % last parameter is always the stage number
    stage=varargin{end};
    v=varargin(1:end-1);

    % penultimate parameter may specify filesystem
    if (~isempty(v) && ischar(v{end}))
        filesystem=v{end};
        v=v(1:end-1);
    else
        filesystem=aap.directory_conventions.remotefilesystem;
    end;

    stagetag=aas_getstagetag(aap,stage);
    doneflagpath=aas_getpath_bydomain(aap,domain,indices,filesystem,stage);
    doneflag=fullfile(doneflagpath, ['done_' stagetag]);

%    aas_cache_put(aap,mfilename,{doneflag doneflagpath stagetag},domain,indices,varargin{:});
end;