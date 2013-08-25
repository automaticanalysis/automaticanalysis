% Automatic analysis
%  Copy files from remote server using rsync+ssh
%  Requires ssh keys to be set up in advance, so that no passwords are
%  required. I typically do this with a host alias in the file ~/.ssh/config
%  as described in this post:
%  http://cusacklabcomputing.blogspot.ca/2012/03/quicker-ssh-config-files-and-automatic.html
%
%  

function [aap]=aas_copyfromremote(aap,host,src,dest,allow404)

% Need to set up garbage collection to limit size of this cache
global aaworker
[pth nme ext]=fileparts(aaworker.parmpath);
cachedir=fullfile(pth,'remotecache');

if (~exist('allow404','var'))
    allow404=false;
end;


allowremotecache=false;
% Only cache single files, if allowed
src=strtrim(src);
if ~any(src==' ')
    allowremotecache=aap.directory_conventions.allowremotecache;
end;

% An aside...
% I learnt something odd about the command "dir"
% If you do exist(pth, 'file') it will find either directories (returns 7) 
% or files (returns a number not 7 or 0 that signifies type). If you do
% exist(pth,'dir') it will only return non-zero if it finds a directory

% Check cache. Either has files, or a flag if nothing came through.
cachehit=false;
if (allowremotecache)    
    % Check to see if already in cache
    md=java.security.MessageDigest.getInstance('MD5');
    md.update(uint8(host),0,length(host));
    md.update(uint8(src),0,length(src));
    md5=md.digest;
    md5(md5<0)=255+md5(md5<0);
    md5=dec2hex(md5);
    md5=md5(:)';
    cachefn=fullfile(cachedir,md5);
    if exist(cachedir,'dir')
        docopy=true;
        if exist(cachefn,'file')
            fprintf('Cache status check\n source: %s\n cache: %s\n dest: %s\n',src,cachefn,dest);
            % Check whether it is in the destination location already
            if exist(dest,'dir')
                [pth nme ext]=fileparts(src);
                dest=fullfile(dest,[nme ext]);
                fprintf(' dest is directory, now:%s\n',dest);
            end;
            if exist(dest,'file') 
                fprintf(' check whether recopy necessary\n');
                stats_src=dir(cachefn)
                stats_dest=dir(dest)
                if (stats_src.bytes==stats_dest.bytes) && (stats_src.datenum==stats_dest.datenum)
                    fprintf(' no copy necessary\n');
                    docopy=false;
                end;
            end;
            % Do the copying if necessary
            if (docopy)
                fprintf(' doing copy\n');
                copyfile(cachefn,dest);
            end;
            cachehit=true;
            aas_log(aap,false,sprintf('Cache of remote file retrieve - got %s',cachefn));
        end;
        if exist([cachefn '_noexist'],'file')
            aas_log(aap,false,sprintf('Cache of remote file retrieve - did not exist\n src:%s\n cache:%s\n dest: %s',src,cachefn,dest));
            cachehit=true;
        end;
    end;
end;

if (~cachehit)
    % Exponential back off if connection refused
    retrydelay=[1 2 4 8 16 32 64 128 256 512 768 1024 2048 1];
    for retry=retrydelay
        % -t option preserves timestamp of remote file
        cmd=sprintf('rsync -t %s:''%s'' %s',host,src,dest);
        [s w]=aas_shell(cmd,allow404);
        if (s==0)
            aas_log(aap,false,sprintf('Retrieved %s from %s',src,host));
            break;
        end;
        if (allow404 && ~isempty(strfind(w,'No such file or directory')))
            break;
        end;
        if (~isempty(strfind(w,'Connection refused')))
            aas_log(aap,false,sprintf('Connection refused in transfer, retry in %d s',retry));
            pause(retry)
        else
            break;
        end;
    end;
    
    % Copy to cache
    if allowremotecache
        aas_makedir(aap,cachedir)
        if exist(dest,'dir')
            [pth nme ext]=fileparts(src);
            dest_exact=fullfile(dest,[nme ext]);
        else
            dest_exact=dest;
        end;
        if (~exist(dest_exact,'file'))
            fid=fopen([cachefn '_noexist'],'w');
            fprintf(fid,sprintf('%f',now));
            fclose(fid);
        else
            copyfile(dest_exact,cachefn);
        end;
    end;
    
    
end;

