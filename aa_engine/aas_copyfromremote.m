function [aap]=aas_copyfromremote(aap,host,src,dest,allow404)
 
if (~exist('allow404','var'))
    allow404=false;
end;


allowremotecache=false;
% Only cache single files, if allowed
src=strtrim(src);
if ~any(src==' ')
    allowremotecache=aap.directory_conventions.allowremotecache;
end;

% Check cache. Either has files, or a flag if nothing came through.
cachehit=false;
if (allowremotecache)
    cachedir='/tmp/remotecache';
    % Check to see if already in cache
    md=java.security.MessageDigest.getInstance('MD5');
    md.update(uint8(host),0,length(host));
    md.update(uint8(src),0,length(src));
    md5=md.digest;
    md5(md5<0)=255+md5(md5<0);
    md5=dec2hex(md5);
    md5=md5(:)';
    cachefn=fullfile(cachedir,md5);
    if (exist(cachedir,'dir'))
        if (exist(cachefn,'file'))
            copyfile(cachefn,dest);
            cachehit=true;
            aas_log(aap,false,sprintf('Cache of remote file retrieve - got %s',cachefn));
        end;
        if exist([cachefn '_noexist'],'file')
            aas_log(aap,false,sprintf('Cache of remote file retrieve - did not exist %s',cachefn));
            cachehit=true;
        end;
    end;
end;

if (~cachehit)
    % Exponential back off if connection refused
    retrydelay=[1 2 4 8 16 32 64 128 256 512];
    for retry=retrydelay
        cmd=sprintf('rsync %s:''%s'' %s',host,src,dest);
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
        if (dest(end)=='/')
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

