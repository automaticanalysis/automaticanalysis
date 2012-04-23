% aa - get path to done flag
% examples:
%  pth=aas_doneflag_getpath(aap,4);  % done flag for stage 4
%  pth=aas_doneflag_getpath(aap,'s3',4);   % path on s3 
%  pth=aas_doneflag_getpath(aap,1,2,4);      % subj 1 sess 2 path for module 4
%  pth=aas_doneflag_getpath(aap,2,4);      % subj 2 path for module 4
% If not specified, remote filesystem defaults to aap.directory_conventions.remotefilesystem
% This differs from the behaviour of aas_getstudypath and derivations,
%  which default to local filesystem

function [doneflag doneflagpath stagetag]=aas_doneflag_getpath(aap,varargin)
global aaworker

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
doneflagname=['done_' stagetag];


% how many parameters left?
switch (length(v))
    case 0
        doneflagpath=aas_getstudypath(aap,filesystem,stage);
    case 1
        doneflagpath=aas_getsubjpath(aap,v{1},filesystem,stage);
    case 2
        doneflagpath=aas_getsesspath(aap,v{1},v{2},filesystem,stage);
end;

doneflag=fullfile(doneflagpath,doneflagname);

% fprintf('Done flag %s\n',doneflag);