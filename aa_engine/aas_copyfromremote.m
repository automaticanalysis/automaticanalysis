% Automatic analysis
%  Copy files from remote server using rsync+ssh
%  Requires ssh keys to be set up in advance, so that no passwords are
%  required. I typically do this with a host alias in the file ~/.ssh/config
%  as described in this post:
%  http://cusacklabcomputing.blogspot.ca/2012/03/quicker-ssh-config-files-and-automatic.html
%
%  Added 'allowcache' input to permit caching on a per-call basis.  Useful if
%  we are adding streams from multiple remote AA locations: some that
%  we want to cache, and some that we don't (e.g., AA analyses on the local
%  machine). CW - 2014/03/05
%
%

function [aap]=aas_copyfromremote(aap,host,src,dest,varargin)

% Set the defaults arguments
vdefaults = { ...
    'allow404',     0, [0 1], ...     % 0 = crash on 404s? 1 = allow them?
    'allowcache',   aap.directory_conventions.allowremotecache, [-1 0 1], ...  % -1 = default to aap.directory_conventions.allowremotecache; 0 = force no; 1 = force yes
    'verbose',      logical(aap.options.verbose), [0 1], ...     % Display all those annoying "Retrieved..." messages?
};

vargs = vargParser(varargin, vdefaults);

allow404 = vargs.allow404;

% Only cache single files, if allowed
allowremotecache = false;
src=strtrim(src);
if ~any(src==' ')
    if vargs.allowcache < 0
        allowremotecache = aap.directory_conventions.allowremotecache;
    else
        allowremotecache = vargs.allowcache;
    end
end;

% An aside...
% I learnt something odd about the command "dir"
% If you do exist(pth, 'file') it will find either directories (returns 7)
% or files (returns a number not 7 or 0 that signifies type). If you do
% exist(pth,'dir') it will only return non-zero if it finds a directory

% Check cache. Either has files, or a flag if nothing came through.
cachehit=false;
if allowremotecache
    % Check to see if already in cache
    md=java.security.MessageDigest.getInstance('MD5');
    if ~isempty(host)
        md.update(uint8(host),0,length(host));
    end
    md.update(uint8(src),0,length(src));
    md5=md.digest;
    md5(md5<0)=255+md5(md5<0);
    md5=dec2hex(md5);
    md5=md5(:)';
    
    % Need to set up garbage collection to limit size of this cache
    global aaworker
    cachedir=fullfile(fileparts(aaworker.parmpath),'remotecache');
    
    cachefn=fullfile(cachedir,md5);
    if exist(cachedir,'dir')
        docopy=true;
        if exist(cachefn,'file')
            aas_log(aap,false,sprintf('Cache status check\n source: %s\n cache: %s\n dest: %s',src,cachefn,dest));
            % Check whether it is in the destination location already
            if exist(dest,'dir')
                [junk, nme, ext]=fileparts(src);
                dest=fullfile(dest,[nme ext]);
                aas_log(aap,false,sprintf(' dest is directory, now:%s',dest));
            end;
            if exist(dest,'file')
                aas_log(aap,false,sprintf(' check whether recopy necessary'));
                stats_src=dir(cachefn);
                stats_dest=dir(dest);
                if (stats_src.bytes==stats_dest.bytes) && (stats_src.datenum==stats_dest.datenum)
                    aas_log(aap,false,sprintf(' no copy necessary'));
                    docopy=false;
                end;
            end;
            % Do the copying if necessary
            if (docopy)
                aas_log(aap,false,sprintf(' doing copy'));
                copyfile(cachefn,dest);
            end;
            aas_log(aap,false,sprintf('Cache of remote file retrieve - got %s',cachefn));
            cachehit=true;
        end;
        if exist([cachefn '_noexist'],'file')
            aas_log(aap,false,sprintf('Cache of remote file retrieve - did not exist\n src:%s\n cache:%s\n dest: %s',src,cachefn,dest));
            cachehit=true;
        end
    end
end

if ~cachehit
    % Exponential back off if connection refused
    retrydelay=[1 2 4 8 16 32 64 128 256 512 768 1024 2048 1];
    for retry=retrydelay
        if isempty(host)
            copyfile(src, dest);
            if vargs.verbose
                 aas_log(aap,false,sprintf('Retrieved %s from %s',src,host),'m');
           end
       else
            % -t option preserves timestamp of remote file
            cmd = sprintf('rsync -t %s:''%s'' %s',host,src,dest);
            [s, w]=aas_shell(cmd,~vargs.verbose,~allow404);
            if (s==0)
                if vargs.verbose
                    aas_log(aap,false,sprintf('Retrieved %s from %s',src,host),'m');
                end
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
        end     
    end;
    
    % Copy to cache
    if allowremotecache
        aas_makedir(aap,cachedir);
        if exist(dest,'dir')
            [junk, nme, ext]=fileparts(src);
            dest=fullfile(dest,[nme ext]);
        end
        if ~exist(dest,'file')
            fid=fopen([cachefn '_noexist'],'w');
            fprintf(fid,sprintf('%f',now));
            fclose(fid);
        else
            copyfile(dest,cachefn);
        end;
    end;
    
    
end;

% Check whether it is an aap_parameter and conversion is required
if strfind(dest,'aap_parameters')
    aap_dest = load(dest);
    if ~isfield(aap_dest.aap.acq_details.subjects(1),'subjname')
        aap_dest.aap = aa_convert_subjects(aap_dest.aap);
        save(dest,'-struct','aap_dest');
        aas_log(aap,false,sprintf('Remote aap %s has been converted',dest),'m');
    end
end
