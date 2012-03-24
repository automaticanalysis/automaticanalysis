% function [aap md5_base64 md5_raw]=aas_md5(aap,filelist[,rootdirectory])
% optional parameter rootdirectory is path added to front of each file
% optional parameter checksumtype = 'data' | 'filestats'
function [aap md5_base64 md5_raw]=aas_md5(aap,filelist,rootdirectory,checksumtype)

if (~exist('checksumtype','var'))
    checksumtype='data';
end;

if (~iscell(filelist))
    tmp=[];
    for ind=1:size(filelist,1);
        tmp{ind}=deblank(filelist(ind,:));
    end;
    filelist=tmp;
end;

% Loop across files
md=java.security.MessageDigest.getInstance('MD5');
for ind=1:length(filelist)
    if (exist('rootdirectory','var') && ~isempty(rootdirectory))
        pth=fullfile(rootdirectory,filelist{ind});
    else
        pth=filelist{ind};
    end;
    
    switch (checksumtype)
        case 'data'
            fid=fopen(pth,'r');

            if (fid==-1)
                aas_log(aap,true,sprintf('Cannot open input file %s',pth));
            end;
            while(~feof(fid))
                data=fread(fid,1024*100,'uint8');
                % [Conor Wild & AVG fix]
                if ~isempty(data)
                    md.update(uint8(data),0,length(data))
                end
            end;
            fclose(fid);
        case 'filestats'
            filestat=dir(pth);
            if (length(filestat)==1)
                data=[num2str(filestat.date) num2str(filestat.bytes)];
                md.update(uint8(data),0,length(data))
            end;
    end;
end;
md5_raw=md.digest;
b64_encoder=sun.misc.BASE64Encoder();
md5_base64=deblank(char(b64_encoder.encode(md5_raw)));